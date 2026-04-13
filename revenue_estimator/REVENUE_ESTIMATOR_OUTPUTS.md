# Revenue Estimator — Output Attribute Reference

This document describes every field in the four output tables written to `revenue_estimation.gpkg` and the companion `revenue_route_direction_summary.csv`.

---

## How the tool works (summary)

The estimator reads raw onboard survey instances and their recorded stop events (boardings and alightings at each stop). For each surveyed trip it:

1. Orders stops by timestamp and projects each stop's GPS point onto the trip geometry to get a distance fraction (0–1 along the route).
2. Reconstructs the passenger load profile — how many people were onboard between each pair of stops.
3. Infers OD flow blocks: groups of passengers assigned a common origin and destination, using the selected behaviour model (FIFO, proportional, or distance-weighted).
4. Assigns a fare to each flow block using a banded distance-fraction model derived from the full terminal-to-terminal fare.
5. Sums revenue per surveyed trip, then scales up to interval-level and daily estimates using headway data from `transit_trips_intervals`.

---

## Table 1: `revenue_stop_profile`

One row per raw stop event per surveyed trip instance. This is the enriched version of `raw.stops` — the original enumerator observations with distance and load fields added.

| Field | Type | Description |
|---|---|---|
| `gid` | integer | Raw stop identifier from RouteLab |
| `name` | text | Stop name as recorded by the enumerator |
| `observer_id` | text | RouteLab observer ID for this stop record |
| `onboard_instance_observer_id` | text | Links this stop to its parent survey instance |
| `stop_time` | timestamp | When the enumerator recorded this stop (renamed from `created_at`). Used to order stops within a trip |
| `board` | float | Passengers boarding at this stop |
| `alight` | float | Passengers alighting at this stop |
| `board_male` / `board_female` | float | Gender-disaggregated boarding counts where recorded |
| `alight_male` / `alight_female` | float | Gender-disaggregated alighting counts where recorded |
| `stop_sequence` | integer | Order of this stop within the trip, derived from `stop_time` ascending |
| `distance_fraction` | float | Where this stop falls along the trip geometry as a normalised fraction (0 = route start, 1 = route end). Computed by projecting the stop GPS point onto the trip linestring |
| `distance_along_m` | float | Same position expressed in metres: `distance_fraction × trip_len_m` |
| `snap_distance_native` | float | Distance in degrees (EPSG:4326) between the raw stop GPS point and its nearest point on the trip line. A large value indicates the stop was recorded far off-route — useful as a QA signal |
| `load_before` | float | Passengers onboard *before* this stop's activity. Equal to `load_after` of the previous stop |
| `load_after` | float | Passengers onboard *after* boarding and alighting: `load_before + board − alight` |
| `trip_ref` | text | Observer trip ID this instance was surveying |
| `route_id` | text | Route identifier from `transit_trips_view` |
| `direction_id` | integer | Direction (0 or 1) from `transit_trips_view` |

---

## Table 2: `revenue_od_matrix`

One row per inferred passenger flow block. Rather than individual passengers, each row represents a cohort of passengers estimated to have travelled from a given origin stop to a given destination stop within a single survey instance.

| Field | Type | Description |
|---|---|---|
| `onboard_instance_id` | text | Survey instance this flow was inferred from |
| `trip_ref` | text | Observer trip ID |
| `trip_gid` | integer | Internal trip ID (links to `transit_trips_view.gid`) |
| `route_id` | text | Route identifier |
| `direction_id` | integer | Direction (0 or 1) |
| `origin_stop_seq` | integer | Stop sequence number where this cohort boarded |
| `destination_stop_seq` | integer | Stop sequence number where this cohort alighted. If `inferred_terminal_alight = True`, this is set to `last_observed_stop + 1`, representing the route terminal |
| `origin_distance_m` | float | Distance along the route at the boarding stop |
| `destination_distance_m` | float | Distance along the route at the alighting stop |
| `traveled_m` | float | Distance this cohort travelled: `destination_distance_m − origin_distance_m` |
| `traveled_fraction` | float | Share of the full route this cohort rode: `destination distance_fraction − origin distance_fraction` (clamped to 0–1) |
| `passenger_count` | float | Number of passengers in this flow block. **See note on fractional values below** |
| `fare_band` | text | Which distance band this journey falls into: `<= 0.33`, `<= 0.66`, or `> upper_band`. Thresholds are configurable at runtime |
| `fare_share` | float | Fraction of the full fare paid by this cohort: 0.50 (band 1), 0.75 (band 2), or 1.00 (band 3) by default |
| `full_fare` | float | Full terminal-to-terminal fare for this trip |
| `estimated_paid_fare` | float | Fare estimated to have been paid per person: `full_fare × fare_share` |
| `flow_revenue` | float | Revenue attributed to this flow block: `estimated_paid_fare × passenger_count` |
| `behavior_model` | text | Allocation model used: `fifo`, `proportional`, or `distance_weighted` |
| `inferred_terminal_alight` | boolean | `True` if these passengers were not observed alighting mid-route but were carried to the terminal. Their alighting was inferred at trip end to recover revenue that would otherwise be lost when enumerators disembark before the final stop |

### Note on fractional `passenger_count`

When the proportional or distance-weighted behaviour model is used, alightings at each stop are distributed across all onboard cohorts in proportion to their size (or distance travelled). This produces mathematically non-integer counts — for example, 2 alightings split across cohorts of 5, 1, and 1 yields 1.429, 0.286, and 0.286.

Passengers cannot logically be fractioned. The fractions are an artefact of the allocation method. However, `flow_revenue` is still internally consistent — summing `flow_revenue` across all rows for a given instance correctly equals `observed_trip_revenue` in the trip trace. The fractions only become a problem if `passenger_count` is used directly as a demand count rather than as a revenue weight. If that is the intended use, rounding at allocation time (using largest-remainder method to preserve totals) should be applied.

---

## Table 3: `revenue_trip_trace`

One row per surveyed trip instance. The per-instance revenue summary, enriched with interval and headway data for scaling.

| Field | Type | Description |
|---|---|---|
| `onboard_instance_id` | text | The specific survey instance |
| `trip_ref` | text | Observer trip ID |
| `trip_gid` | integer | Internal trip ID |
| `route_id` | text | Route identifier |
| `direction_id` | integer | Direction (0 or 1) |
| `vehicle_name` | text | Operator or vehicle type (e.g. Taxi Collectif, STS bus) |
| `origin` | text | Origin terminal name |
| `destination` | text | Destination terminal name |
| `interval_observer_id` | text | Raw interval FK from the onboard instance (Postgres path). Null in GeoPackage mode — interval is resolved via `interval_name` instead |
| `interval_name` | text | Time-of-day interval the survey was conducted in (e.g. `morning_peak`, `afternoon`) |
| `full_fare` | float | Full terminal-to-terminal fare for this trip |
| `total_board` | float | Sum of all boardings recorded across all stops in this instance |
| `total_alight` | float | Sum of all alightings recorded across all stops. Does not include terminal alightings — those are recorded separately in `terminal_alight` |
| `terminal_alight` | float | Passengers carried to the terminal without a recorded alighting, recovered by the terminal inference step. These contribute to `observed_trip_revenue` but are not in `total_alight` |
| `final_load` | float | Passengers still onboard after the last recorded stop: `total_board − total_alight`. A non-zero value means the enumerator disembarked before the terminal. After the terminal inference fix, these passengers are captured in `terminal_alight` and their revenue included in `observed_trip_revenue` |
| `observed_trip_revenue` | float | Total revenue attributed to this single survey run, summed across all OD flow blocks including terminal alightings |
| `avg_inferred_trip_fraction` | float | Weighted average of `traveled_fraction` across all OD rows for this instance, weighted by `passenger_count`. Indicates on average how much of the route passengers rode (0 = none, 1 = full route) |
| `avg_inferred_fare` | float | Weighted average fare paid per passenger. Compare to `full_fare` to see what fraction of the full fare the average passenger actually paid |
| `behavior_model` | text | Allocation model used |
| `fare_model_type` | text | Always `banded_fraction_of_full_fare` — documents the fare model in use |
| `interval_gid` | integer | Internal interval ID joined from `transit_intervals` |
| `start_time` / `end_time` | time | Start and end time of the interval |
| `active` | integer | Whether this interval is active (1) or not (0) |
| `headway_secs` | float | Average headway in seconds for this trip in this interval, from `transit_trips_intervals` |
| `headway_estimation_method` | text | How the headway was derived: `from_own_freq_surveys` (direct observation), `from_similar_agency_interval` (ratio-based imputation from another interval of the same operator), or `from_user_default` (manual fallback value) |
| `interval_duration_secs` | float | Length of the interval in seconds (`end_time − start_time`) |
| `vehicle_trips_in_interval` | float | Estimated number of vehicle trips of this type during the interval: `interval_duration_secs ÷ headway_secs` |
| `estimated_interval_revenue` | float | Total estimated revenue for all trips of this type in this interval, assuming each trip earns the same as this surveyed instance: `observed_trip_revenue × vehicle_trips_in_interval` |

---

## Table 4: `revenue_route_direction_summary`

One row per route × direction combination. The final aggregated output for planning and reporting use. Also written as `revenue_route_direction_summary.csv`.

| Field | Type | Description |
|---|---|---|
| `route_id` | text | Route identifier |
| `direction_id` | integer | Direction (0 or 1) |
| `origin` | text | Origin terminal name |
| `destination` | text | Destination terminal name |
| `vehicle_name` | text | Operator or vehicle type |
| `full_fare` | float | Full terminal-to-terminal fare |
| `sample_trip_count` | integer | Number of distinct surveyed instances for this route/direction |
| `sample_total_board` | float | Total boardings across all survey instances |
| `sample_avg_trip_revenue` | float | Average `observed_trip_revenue` across all surveyed instances, regardless of interval. A rough overall benchmark |
| `headway_secs__<interval>` | float | Average headway for this route/direction in the named interval |
| `vehicle_trips__<interval>` | float | Estimated number of vehicle trips in the named interval: `interval_duration ÷ headway_secs` |
| `sample_avg_trip_revenue__<interval>` | float | Average observed revenue from surveys conducted specifically in the named interval. This is the key input to interval-level scaling. **Null if no surveys were conducted in this interval** |
| `interval_revenue__<interval>` | float | Estimated total revenue for this route/direction during the named interval: `sample_avg_trip_revenue__<interval> × vehicle_trips__<interval>`. Null if `sample_avg_trip_revenue__<interval>` is null |
| `active__<interval>` | float | Whether this interval is active (1) or not (0) |
| `daily_estimated_revenue` | float | Sum of `interval_revenue__*` across all intervals — the total estimated daily revenue for this route/direction |
| `active_interval_count` | float | Number of active intervals that contributed a non-null revenue estimate to the daily total |

### Important caveat on `daily_estimated_revenue`

`daily_estimated_revenue` is the sum of interval revenues. Any interval with no survey coverage contributes nothing to this total — its `sample_avg_trip_revenue__<interval>` is null and its `interval_revenue__<interval>` is excluded from the sum. If a route was only surveyed during morning peak but also operates in the afternoon, the afternoon interval contributes zero to the daily estimate even though `vehicle_trips__afternoon` is populated. This is by design — the tool will not fabricate revenue for unobserved intervals — but means `daily_estimated_revenue` should be interpreted as a lower bound when survey coverage is incomplete.

---

## QA tables

### `revenue_qa_trips`

One row per QA issue raised during processing. Issues are raised at the instance level.

| Field | Description |
|---|---|
| `onboard_instance_id` | Instance where the issue was detected |
| `trip_ref` | Associated trip |
| `route_id`, `direction_id` | Route context |
| `qa_issue_code` | Machine-readable issue code (see below) |
| `qa_issue_detail` | Human-readable description with specific values |
| `severity` | `warning` (processing continued) or `error` (data quality issue) |

**Issue codes:**

| Code | Meaning |
|---|---|
| `NO_STOPS` | No raw stop events found for this onboard instance |
| `MISSING_TRIP` | The instance's `trip_ref` was not found in `transit_trips_view` |
| `NEGATIVE_INTERMEDIATE_LOAD` | Alightings exceeded the onboard count at some point — more people got off than were on. Usually caused by recording errors |
| `NONZERO_FINAL_LOAD` | Passengers remained onboard after the last recorded stop (non-zero `final_load`). After the terminal inference fix, these passengers are recovered into `terminal_alight` rather than lost |
| `DUPLICATE_STOP_TIMESTAMPS` | Two or more stops share the same timestamp — stop ordering may be unreliable |
| `ALIGHTING_SHORTFALL` | The allocation model could not assign all alightings at a stop to existing cohorts — more people alighted than were onboard |
| `UNMATCHED_ONBOARD` | After terminal inference, passengers remained unallocated. Should not occur under normal conditions |
| `PROCESSING_FAILED` | An unhandled exception occurred during processing of this instance. The instance was skipped |

### `revenue_qa_summary`

Aggregate summary of QA issues across all instances.

| Field | Description |
|---|---|
| `qa_issue_code` | Issue code |
| `trip_count` | Number of distinct instances affected |
| `record_count` | Total number of QA records with this code |

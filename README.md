# TfC Tools QGIS Plugin

**TfC Tools** is a QGIS plugin suite developed by [Transport for Cairo (TfC)](https://transportforcairo.com) to standardize and streamline transport data processing workflows.  
It provides tools commonly used in TfC’s research and transport data management projects, covering field data migration, GTFS generation, and flow estimation.

---

## Plugin Overview

TfC Tools includes three Processing Toolbox plugins that can be used independently or as part of a workflow:

**Workflow:**  
`RouteLab → RL2SDI / RL2GPKG / SDI2GPKG → GeoPackage / SDI → GIS2GTFS / Revenue Estimator / Vehicle & Passenger Flow`

## Plugin Structure

TfC Tools v2.0 introduces two tool groups:

### 1. RouteLab Tools
Tools that migrate field survey data from TfC’s RouteLab database into standardized tables, either in GeoPackage format or in a PostgreSQL/PostGIS Spatial Data Infrastructure (SDI).

| Tool | Purpose |
|-----|-----|
| Export RouteLab to GeoPackage | Exports RouteLab survey data into a portable GeoPackage format for analysis without requiring a database. |
| Export SDI (PostGIS) to GeoPackage | Exports an existing SDI database to a GeoPackage format. |
| RL2SDI (RouteLab to SDI Migration) | Migrates RouteLab survey data into a standardized PostGIS SDI schema. |

### 2. GIS Tools
Tools for transport network analysis, GTFS generation, and transport metrics.

| Tool | Purpose |
|-----|-----|
| GIS2GTFS | Converts GIS data in a standards schema, either from a PostGIS SDI or a Geopackage file, into a valid GTFS feed. |
| Refresh SDI derived layers | Rebuilds derived tables and materialized views in an SDI or GeoPackage. |
| Trip and Route Revenue Estimator | Estimates potential revenue for routes and trips based on ridership and fare inputs. |
| Vehicle and Passenger Flow Estimator | Estimates vehicle flows and passenger flows on road segments using standardized SDI or GeoPackage data. |
---

## Installation

### Option 1 — Manual (from ZIP)
1. Download the latest plugin ZIP from this GitHub repository.  
2. In QGIS, go to **Plugins → Manage and Install Plugins…**  
3. Select **Install from ZIP**, browse to the downloaded file, and click **Install Plugin**.

### Option 2 — From QGIS Plugin Repository
Once approved, TfC Tools will also be available directly through QGIS:
1. Open **Plugins → Manage and Install Plugins…**  
2. Search for “TfC Tools”.  
3. Click **Install Plugin**.

After installation, ensure **Plugins → TfC Tools** is checked, and access the tools from the **Processing Toolbox** panel.

---

## 1. RL2SDI Migration Plugin

**Purpose:**  
Migrates field survey data from TfC’s RouteLab database into a standardized SDI (PostGIS).

**Inputs:**
- RouteLab database connection (credentials provided by TfC).  
- PostGIS SDI connection for data migration.  
- Project ID (unique identifier per RouteLab project).
- Headways (optional, for trips with missing headway data).

**Output:**  
Creates two new schemas in the target database:
- `raw` — cleaned field data from RouteLab  
- `transit` — processed, structured results

A complete list of output tables and fields for both schemas is included in the User Guide.

**Next Steps:**  
Output can be analyzed directly in QGIS or used as input to the *GIS2GTFS* or *Vehicle and Passenger Flow* plugins.

---

## 2. Export RouteLab to GeoPackage

**Purpose:**  
Exports RouteLab field survey data into a GeoPackage file following TfC’s standardized schema.  
This allows working with the data without requiring a PostgreSQL database.

**Inputs:**
- RouteLab database connection (credentials provided by TfC).  
- Project ID (unique identifier per RouteLab project).  
- Headways (optional, for trips with missing headway data).
- Output folder.

**Outputs:**
- A GeoPackage file containing standardized datasets equivalent to the `raw` and `transit` structures.

**Use Cases:**
- Offline analysis without database access  
- Data sharing across teams  
- Input to other TfC Tools (GIS2GTFS, Flow, Revenue)

---

## 3. Export SDI (PostGIS) to GeoPackage

**Purpose:**  
Exports an existing PostgreSQL/PostGIS SDI database into a GeoPackage format.

**Inputs:**
- PostgreSQL SDI connection (standard schema).  
- Headways (optional, for trips with missing headway data).
- Output folder.

**Outputs:**
- A GeoPackage replicating the SDI schema and datasets.

**Use Cases:**
- Sharing SDI datasets without database access  
- Creating portable project files  
- Backup and archiving

---

## 4. GIS2GTFS Plugin

**Purpose:**  
Generates a GTFS feed from a standardized dataset stored in either a PostGIS SDI or a GeoPackage.

**Inputs:**
- SDI connection or GeoPackage (standard schema from *RL2SDI* or *RL2GPKG* or following the schema described in the User Guide). 
- Feed version, start and end dates, service ID.  
- Option to enable/disable continuous drop-off/pick-up.  
- Two output folders:  
  - **Raw folder:** temporary files (for internal use)  
  - **Output folder:** final GTFS `.txt` files

**Outputs:**
- A complete GTFS dataset (9 text files) exported to the selected folder.

**Validation (recommended):**
Validate using [MobilityData GTFS Validator](https://github.com/MobilityData/gtfs-validator):
1. Zip the `.txt` files.  
2. Run the validator and review results.  
3. Correct any *errors* before distribution; *warnings* are optional.

---

## 5. Refresh SDI derived layers

**Purpose:**  
Rebuilds derived tables and materialized views in the TfC SDI schema.

**Inputs:**
- SDI database or GeoPackage (standard schema).  
- Optional setting to include QA layers.

**Outputs:**
- Updated derived tables (e.g. trips_view, trip_stops_sequence, od_stats).

**Notes:**
- Existing layers are overwritten.  
- Missing optional layers are skipped with warnings.

---

## 6. Trip and Route Revenue Estimator

**Purpose:**  
Estimates revenue for transport routes and trips based on observed demand and fare assumptions.

**Inputs:**
- SDI database or GeoPackage (standard schema from *RL2SDI* or *RL2GPKG* or following the schema described in the User Guide).  
- Fare assumptions.  
- Output folder.

**Outputs:**
- Tables summarizing estimated revenue per trip and per route.

**Use Cases:**
- Service performance analysis  
- Scenario testing for fares  
- Planning and financial assessment

---

## 7. Vehicle and Passenger Flow Plugin

**Purpose:**  
Estimates vehicle and passenger flows per road segment using standardized SDI or GeoPackage data.

**Inputs:**
- PostGIS connection or GeoPackage (standard schema from *RL2SDI* or *RL2GPKG* or following the schema described in the User Guide).  
- Output folder for generated files.

**Process Summary:**
1. Builds trip segments from stop sequences.  
2. Uses onboard survey data and OD statistics.  
3. Applies headway information to estimate vehicle supply. 
4. Calculates vehicle and passenger volumes per segment and time interval.  
5. Outputs the results for visualization in QGIS.

**Outputs:**
- `vehicle_passenger_flow.gpkg` — Vehicle and passenger flow layers  

---

## Technical Notes

- The plugins rely on PostGIS connections configured in QGIS.  
- Compatible with **QGIS 3.40+**.
- Designed for reproducible workflows using TfC’s standard SDI schema.  

---

## Documentation

Full details, figures, and database schema descriptions are available in the **User Guide (PDF)**:  
[**tfc_tools_user_guide.pdf**](https://github.com/transportforcairo/tfc_tools/blob/main/tfc_tools_user_guide.pdf)


---

## License

TfC Tools is released under the **GNU General Public License v2 (GPL-2.0)**.  
See the [LICENSE](./LICENSE) file for details.

---

## Version

**Current version:** 2.0  |  **Release date:** April 2026

---

## Contact

**Transport for Cairo (TfC)**  
[transportforcairo.com](https://transportforcairo.com)  
For technical questions or contributions, please open an [issue](https://github.com/transportforcairo/tfc_tools/issues).

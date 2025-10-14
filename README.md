# TfC Tools QGIS Plugin

**TfC Tools** is a QGIS plugin suite developed by [Transport for Cairo (TfC)](https://transportforcairo.com) to standardize and streamline transport data processing workflows.  
It provides tools commonly used in TfC’s research and transport data management projects, covering field data migration, GTFS generation, and flow estimation.

---

## Plugin Overview

TfC Tools includes three Processing Toolbox plugins that can be used independently or as part of a workflow:

**Workflow:**  
`RouteLab → RL2SDI → GIS2GTFS → Vehicle and Passenger Flow`

| Plugin | Purpose |
|---------|----------|
| **RL2SDI (RouteLab to SDI Migration)** | Migrates RouteLab field survey data into a standardized PostGIS Spatial Data Infrastructure (SDI). |
| **GIS2GTFS** | Converts a PostGIS SDI (in TfC’s standard schema) into a valid GTFS feed. |
| **Vehicle and Passenger Flow** | Estimates vehicle and passenger flows on road segments based on GTFS data. |

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

**Output:**  
Creates two new schemas in the target database:
- `raw` — cleaned field data from RouteLab  
- `transit` — processed, structured results

**Next Steps:**  
Output can be analyzed directly in QGIS or used as input to the *GIS2GTFS* or *Vehicle and Passenger Flow* plugins.

---

## 2. GIS2GTFS Plugin

**Purpose:**  
Generates a GTFS feed from a standardized SDI in PostGIS.

**Inputs:**
- SDI connection (from *RL2SDI*).  
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

## 3. Vehicle and Passenger Flow Plugin

**Purpose:**  
Estimates vehicle and passenger flows per road segment based on GTFS data.

**Inputs:**
- GTFS `.zip` file (e.g. output from *GIS2GTFS*).  
- PostGIS connection (standard schema from *RL2SDI*).  
- Output folder for generated files.

**Process Summary:**
1. Reads GTFS and extracts the geographic boundary.  
2. Downloads and filters the OSM road network within that boundary.  
3. Computes Fréchet distance between trip and road segments.  
4. Calculates vehicle and passenger volumes per segment and time interval.  
5. Outputs the results for visualization in QGIS.

**Outputs:**
- `veh_flow.gpkg` — Vehicle flow (morning and afternoon peaks)  
- `passenger_flow.gpkg` — Passenger flow (morning and afternoon peaks)  
- `debug_boundary.geojson` — boundary visualization layer

**Output Fields:**
| Field | Description |
|--------|-------------|
| `fid` | Auto-generated unique ID |
| `gid` | Road segment ID from OSM |
| `interval_name` | Time interval name |
| `value` | Estimated passengers or vehicles per segment per hour |

---

## Technical Notes

- The plugins rely on PostGIS connections configured in QGIS.  
- Compatible with **QGIS 3.0 and later** (tested on 3.34).  
- Designed for reproducible workflows using TfC’s standard SDI schema.  

---

## Documentation

A detailed **User Guide** (with screenshots and schema tables) will be added to this repository soon.

---

## License

TfC Tools is released under the **GNU General Public License v2 (GPL-2.0)**.  
See the [LICENSE](./LICENSE) file for details.

---

## Contact

**Transport for Cairo (TfC)**  
[transportforcairo.com](https://transportforcairo.com)  
For technical questions or contributions, please open an [issue](https://github.com/transportforcairo/tfc_tools/issues).

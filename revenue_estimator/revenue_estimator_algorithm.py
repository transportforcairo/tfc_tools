# -*- coding: utf-8 -*-
from __future__ import annotations

import os

from ..tfc_tools_common import ensure_paths, ensure_deps
ensure_paths()

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterEnum,
    QgsProcessingParameterProviderConnection,
    QgsProcessingParameterFile,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterDefinition,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
)

from ..tfc_tools_common.sdi_io import SDISource
from .engine import RevenueEstimator, RevenueConfig

# for icon
from qgis.PyQt.QtGui import QIcon

def _icon_path(*parts):
    # __file__ is inside tfc_tools/rl2sdi/
    here = os.path.dirname(__file__)
    root = os.path.abspath(os.path.join(here, ".."))     # tfc_tools/
    return os.path.join(root, *parts)


class RevenueEstimatorAlgorithm(QgsProcessingAlgorithm):
    SDI_SOURCE = "SDI_SOURCE"
    SDI_CONNECTION = "SDI_CONNECTION"
    SDI_GPKG = "SDI_GPKG"
    OUTPUT_FOLDER = "OUTPUT_FOLDER"
    BEHAVIOR_MODEL = "BEHAVIOR_MODEL"
    QA_MODE = "QA_MODE"
    DROP_ZERO_STOPS = "DROP_ZERO_STOPS"
    BAND1_MAX = "BAND1_MAX"
    BAND1_SHARE = "BAND1_SHARE"
    BAND2_MAX = "BAND2_MAX"
    BAND2_SHARE = "BAND2_SHARE"
    MIN_HEADWAY_SECS = "MIN_HEADWAY_SECS"
    MAX_SNAP_DISTANCE_M = "MAX_SNAP_DISTANCE_M"
    FINAL_LOAD_TOLERANCE = "FINAL_LOAD_TOLERANCE"
    WRITE_OD_MATRIX = "WRITE_OD_MATRIX"
    WRITE_STOP_PROFILE = "WRITE_STOP_PROFILE"
    WRITE_QA_TABLES = "WRITE_QA_TABLES"

    def tr(self, s: str) -> str:
        return QCoreApplication.translate(self.__class__.__name__, s)

    def name(self):
        return "revenue_estimator"

    def displayName(self):
        return self.tr("Trip and Route Revenue Estimator")

    def group(self):
        return self.tr("02 GIS Tools")

    def groupId(self):
        return "gis_tools"

    def createInstance(self):
        return RevenueEstimatorAlgorithm()

    def shortHelpString(self):
        return """
    <b>Purpose of the Plugin</b>
    This tool estimates revenue for transport services using onboard trip instances, raw stop counts, trip geometry, full route fare, and interval/headway tables.
    The tool works with either:
    • an SDI PostGIS database connection, or
    • an SDI GeoPackage exported from TfC Tools / RouteLab workflows.
    No OD survey input is required. Passenger OD flows are inferred from boarding and alighting counts observed at raw stops.<br>

    <b>What this tool does</b>
    1. Reads onboard trip instances, raw stops, trip reference data, and interval/headway tables.
    2. Sequences stops using raw stop timestamps.
    3. Projects each stop onto the trip geometry to measure cumulative distance traveled along the route.
    4. Reconstructs onboard passenger load along the trip using boarding and alighting counts.
    5. Infers OD flow blocks from aggregate stop counts using the selected behaviour model.
    6. Estimates fare paid by each inferred passenger flow using a distance-based, banded fare model derived from the full route fare.
    7. Calculates observed surveyed-trip revenue, interval-scaled service revenue, and route + direction summary outputs
    8. Writes traceable output tables for QA, stop profiles, OD matrix, trip trace, and summary tables.<br>

    <b>Important assumptions:</b>
    • transit_trips_view.fare is treated as the full terminal-to-terminal fare.
    • Stop sequence is derived from raw_stops.created_at.
    • Distance is measured along trip geometry, not straight-line distance.
    • Revenue is estimated by inferred OD flow blocks, not directly observed passenger-level OD records.<br>

    <b>Outputs</b>
    Outputs include stop profiles, OD matrix, revenue calculation trace, and route + direction summary tables, and optional QA tables.<br>

    <b>Recommended starting setup:</b>
    • Behaviour model = Proportional
    • QA mode = Warn
    • Fare bands:
    - band 1 max fraction = 0.33
    - band 1 fare share = 0.50
    - band 2 max fraction = 0.66
    - band 2 fare share = 0.75<br>

    <b>Documentation</b>
    For full inputs, parameters, outputs, and interpretation notes, refer to the User Guide.
    <a href="https://github.com/transportforcairo/tfc_tools/blob/main/tfc_tools_user_guide.pdf">TfC Tools User Guide</a>
    """

    def icon(self):
        return QIcon(_icon_path("icons", "GIS-icon.svg"))
    # Optional: some QGIS builds prefer this for SVG
    def svgIconPath(self):
        return _icon_path("icons", "GIS-icon.svg")
    

    def initAlgorithm(self, config=None):
        src_param = QgsProcessingParameterEnum(
            self.SDI_SOURCE,
            self.tr("SDI data source"),
            options=["PostGIS (QGIS connection)", "GeoPackage (RouteLab / SDI export)"],
            defaultValue=0,
        )
        src_param.setMetadata({
            "widget_wrapper": {
                "class": "tfc_tools.tfc_tools_common.sdi_source_toggle_wrapper.SDISourceToggleWidgetWrapper",
                "enable_map": {"0": [self.SDI_CONNECTION], "1": [self.SDI_GPKG]},
            }
        })
        self.addParameter(src_param)

        self.addParameter(QgsProcessingParameterProviderConnection(
            self.SDI_CONNECTION,
            self.tr("SDI PostGIS connection (required when SDI data source = PostGIS)"),
            provider="postgres",
        ))
        self.parameterDefinition(self.SDI_CONNECTION).setFlags(
            self.parameterDefinition(self.SDI_CONNECTION).flags() | QgsProcessingParameterDefinition.FlagOptional
        )

        self.addParameter(QgsProcessingParameterFile(
            self.SDI_GPKG,
            self.tr("SDI GeoPackage (required when SDI data source = GeoPackage)"),
            behavior=QgsProcessingParameterFile.File,
            fileFilter="GeoPackage (*.gpkg)",
            optional=True,
        ))

        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_FOLDER, self.tr("Output folder")))

        self.addParameter(QgsProcessingParameterEnum(
            self.BEHAVIOR_MODEL,
            self.tr("Behavior model"),
            options=["FIFO", "Proportional", "Distance-weighted"],
            defaultValue=1,
        ))
        self.addParameter(QgsProcessingParameterEnum(
            self.QA_MODE,
            self.tr("QA mode"),
            options=["Strict", "Warn", "Repair minor issues"],
            defaultValue=1,
        ))
        self.addParameter(QgsProcessingParameterBoolean(self.DROP_ZERO_STOPS, self.tr("Drop zero-board / zero-alight stops"), defaultValue=True))

        self.addParameter(QgsProcessingParameterNumber(self.BAND1_MAX, self.tr("Fare band 1: maximum route fraction"), type=QgsProcessingParameterNumber.Double, defaultValue=0.33, minValue=0.0, maxValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(self.BAND1_SHARE, self.tr("Fare band 1: share of full fare"), type=QgsProcessingParameterNumber.Double, defaultValue=0.50, minValue=0.0, maxValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(self.BAND2_MAX, self.tr("Fare band 2: maximum route fraction"), type=QgsProcessingParameterNumber.Double, defaultValue=0.66, minValue=0.0, maxValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(self.BAND2_SHARE, self.tr("Fare band 2: share of full fare"), type=QgsProcessingParameterNumber.Double, defaultValue=0.75, minValue=0.0, maxValue=1.0))

        for name, label, default in [
            (self.MIN_HEADWAY_SECS, "Minimum valid headway (seconds)", 60),
            (self.MAX_SNAP_DISTANCE_M, "Maximum stop-to-line snap distance (meters, QA only)", 100),
            (self.FINAL_LOAD_TOLERANCE, "Final load tolerance", 0),
        ]:
            p = QgsProcessingParameterNumber(name, self.tr(label), type=QgsProcessingParameterNumber.Double, defaultValue=default, minValue=0.0)
            p.setFlags(p.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
            self.addParameter(p)
        for name, label, default in [
            (self.WRITE_OD_MATRIX, "Write OD matrix", True),
            (self.WRITE_STOP_PROFILE, "Write stop profile", True),
            (self.WRITE_QA_TABLES, "Write QA tables", True),
        ]:
            p = QgsProcessingParameterBoolean(name, self.tr(label), defaultValue=default)
            p.setFlags(p.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
            self.addParameter(p)

    def checkParameterValues(self, parameters, context):
        source_mode = self.parameterAsEnum(parameters, self.SDI_SOURCE, context)
        connection_name = self.parameterAsString(parameters, self.SDI_CONNECTION, context)
        gpkg_path = self.parameterAsFile(parameters, self.SDI_GPKG, context)
        b1m = self.parameterAsDouble(parameters, self.BAND1_MAX, context)
        b2m = self.parameterAsDouble(parameters, self.BAND2_MAX, context)
        b1s = self.parameterAsDouble(parameters, self.BAND1_SHARE, context)
        b2s = self.parameterAsDouble(parameters, self.BAND2_SHARE, context)
        if source_mode == 0:
            if not connection_name:
                return (False, self.tr("Please choose an SDI PostGIS connection."))
        else:
            if not gpkg_path:
                return (False, self.tr("Please select an SDI GeoPackage file."))
            if not os.path.exists(gpkg_path):
                return (False, self.tr("The selected GeoPackage path does not exist."))
            if not gpkg_path.lower().endswith(".gpkg"):
                return (False, self.tr("The selected file is not a .gpkg GeoPackage."))
        if not (0 < b1m < b2m <= 1.0):
            return (False, self.tr("Fare band maximum fractions must satisfy 0 < band 1 max < band 2 max <= 1."))
        if not (0 <= b1s <= b2s <= 1.0):
            return (False, self.tr("Fare shares must satisfy 0 <= band 1 share <= band 2 share <= 1."))
        return super().checkParameterValues(parameters, context)

    def processAlgorithm(self, parameters, context, feedback):
        ensure_deps(show_ui=True)
        source_mode = self.parameterAsEnum(parameters, self.SDI_SOURCE, context)
        if source_mode == 0:
            source = SDISource(mode="postgres", conn_name=self.parameterAsConnectionName(parameters, self.SDI_CONNECTION, context))
        else:
            source = SDISource(mode="gpkg", gpkg_path=self.parameterAsFile(parameters, self.SDI_GPKG, context))

        behavior_map = {0: "fifo", 1: "proportional", 2: "distance_weighted"}
        qa_map = {0: "strict", 1: "warn", 2: "repair_minor"}
        cfg = RevenueConfig(
            behavior_model=behavior_map[self.parameterAsEnum(parameters, self.BEHAVIOR_MODEL, context)],
            qa_mode=qa_map[self.parameterAsEnum(parameters, self.QA_MODE, context)],
            drop_zero_stops=self.parameterAsBool(parameters, self.DROP_ZERO_STOPS, context),
            band1_max_fraction=self.parameterAsDouble(parameters, self.BAND1_MAX, context),
            band1_share=self.parameterAsDouble(parameters, self.BAND1_SHARE, context),
            band2_max_fraction=self.parameterAsDouble(parameters, self.BAND2_MAX, context),
            band2_share=self.parameterAsDouble(parameters, self.BAND2_SHARE, context),
            min_headway_secs=self.parameterAsDouble(parameters, self.MIN_HEADWAY_SECS, context),
            max_snap_distance_m=self.parameterAsDouble(parameters, self.MAX_SNAP_DISTANCE_M, context),
            final_load_tolerance=self.parameterAsDouble(parameters, self.FINAL_LOAD_TOLERANCE, context),
            write_od_matrix=self.parameterAsBool(parameters, self.WRITE_OD_MATRIX, context),
            write_stop_profile=self.parameterAsBool(parameters, self.WRITE_STOP_PROFILE, context),
            write_qa_tables=self.parameterAsBool(parameters, self.WRITE_QA_TABLES, context),
        )
        output_folder = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        engine = RevenueEstimator(source=source, output_folder=output_folder, config=cfg)
        return engine.run(feedback=feedback)

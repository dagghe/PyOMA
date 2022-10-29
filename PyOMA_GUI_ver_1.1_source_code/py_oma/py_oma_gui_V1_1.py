from PyQt5.QtWidgets import QMainWindow, QApplication, QLabel, QPushButton, QFileDialog, QTextBrowser, QMessageBox, \
    QTextEdit, QScrollArea, QComboBox, QFormLayout, QGroupBox, QVBoxLayout, QListWidget, QErrorMessage, QHBoxLayout, \
    QToolButton, QLineEdit, QDialogButtonBox, QDialog, QWidget, QCheckBox, QTableWidget, QTableWidgetItem, QInputDialog
from PyQt5.QtCore import Qt, QDir, QFileInfo
from PyQt5 import uic, QtGui
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import signal
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
import py_oma.PyOMA as oma
import py_oma.drawing_tools_3d as drawing_tools_3d
import sys
import os
import random
import matplotlib
matplotlib.use('QT5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from py_oma.utils import resolve_path
from dataclasses import dataclass, field
import mplcursors
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D


class UI(QMainWindow):


    def __init__(self):
        super(UI, self).__init__()
        
        # Load the ui file + basic properties
        uic.loadUi(resolve_path("py_oma_V1_1.ui"), self)
        self.setWindowIcon(QtGui.QIcon('logopyoma_ico.png'))
        self.pathFolder = None
        self.pathLoadFile = ""
        self.inputData = 0
        self.resultsDirectory = ""
        self.plotWidgetGeometry = FigureCanvas()
        self.layGeometryPlot = QVBoxLayout(self.content_PlotInitialGeometry)
        self.layGeometryPlot.setContentsMargins(0, 0, 0, 0)
        self.layGeometryTool = QVBoxLayout(self.content_ToolbarInitialGeometry)
        self.layGeometryTool.setContentsMargins(0, 0, 0, 0)
        self.plotWidgetFddSvp = FigureCanvas()
        self.layFddSvpPlot = QVBoxLayout(self.content_PlotFddSvp)
        self.layFddSvpPlot.setContentsMargins(0, 0, 0, 0)
        self.layFddSvpTool = QVBoxLayout(self.content_ToolbarFddSvp)
        self.layFddSvpTool.setContentsMargins(0, 0, 0, 0)
        self.plotWidgetFdd = FigureCanvas()
        self.layFddPlot = QVBoxLayout(self.content_PlotFdd)
        self.layFddPlot.setContentsMargins(0, 0, 0, 0)
        self.layFddTool = QVBoxLayout(self.content_ToolbarFdd)
        self.layFddTool.setContentsMargins(0, 0, 0, 0)
        self.plotWidgetFddGeom = FigureCanvas()
        self.layFddPlotGeom = QVBoxLayout(self.content_PlotFddGeom)
        self.layFddPlotGeom.setContentsMargins(0, 0, 0, 0)
        self.layFddToolGeom = QVBoxLayout(self.content_ToolbarFddGeom)
        self.layFddToolGeom.setContentsMargins(0, 0, 0, 0)
        self.plotWidgetSsi = FigureCanvas()
        self.laySsiPlot = QVBoxLayout(self.content_PlotSsi)
        self.laySsiPlot.setContentsMargins(0, 0, 0, 0)
        self.laySsiTool = QVBoxLayout(self.content_ToolbarSsi)
        self.laySsiTool.setContentsMargins(0, 0, 0, 0)
        self.plotWidgetSsiGeom = FigureCanvas()
        self.laySsiPlotGeom = QVBoxLayout(self.content_PlotSsiGeom)
        self.laySsiPlotGeom.setContentsMargins(0, 0, 0, 0)
        self.laySsiToolGeom = QVBoxLayout(self.content_ToolbarSsiGeom)
        self.laySsiToolGeom.setContentsMargins(0, 0, 0, 0)
        self._figuresFDD = {}
        self._figuresSSI = {}
        self.nodesDict = {}
        # self.channelsDict: {'Name of channel': [ID of node, direction (0-->x, 1-->y, 2-->z)]}
        self.channelsDict = {}
        self.channelNamesDict = {}
        self.no_channels = 0
        self.nodes = []
        self.connectivity = []
        self.changed_items = []
        self.changed_items_name = []
        self.parameters = Parameters()
        self.fddhelper = 0
        self.my_clicker_helper = 0
        self.modesDict = {}
        self.setup_clicks = 1
        self.num_setup = None

        # Assign Widgets of Designer to py code - Import Data
        self.labelCreateFolder = self.findChild(QLabel, "label_CreateFolder")
        self.labelNumOfChannels = self.findChild(QLabel, "label_NumChannels")
        self.labelOpenFile = self.findChild(QLabel, "label_OpenFile")
        self.labelTimeSteps = self.findChild(QLabel, "label_TimeSteps")
        self.buttonClearAll = self.findChild(QPushButton, "pushButton_ClearAllSingle")
        self.buttonCreateFolder = self.findChild(QPushButton, "pushButton_CreateFolder")
        self.buttonLoadData = self.findChild(QPushButton, "pushButton_LoadData")
        self.buttonSubmitSetup = self.findChild(QPushButton, "pushButton_SubmitSetup")
        self.tableInputData = self.findChild(QTableWidget, "tableWidget_InputData")

        # Geometry
        self.comboBoxChannelNames = self.findChild(QComboBox, "comboBox_ChannelsNames")
        self.comboBoxChannelAssignments = self.findChild(QComboBox, "comboBox_AllChannelsAssignments")
        self.labelOpenFileConnectivity = self.findChild(QLabel, "label_OpenFileConnectivity")
        self.labelOpenFileNodes = self.findChild(QLabel, "label_OpenFileNodes")
        self.buttonClearChannels = self.findChild(QPushButton, "pushButton_ClearChannels")
        self.buttonCreateGeometry = self.findChild(QPushButton, "pushButton_CreateGeometry")
        self.buttonLoadConnectivity = self.findChild(QPushButton, "pushButton_LoadConnectivity")
        self.buttonLoadNodes = self.findChild(QPushButton, "pushButton_LoadNodes")
        self.buttonSubmitChannels = self.findChild(QPushButton, "pushButton_SubmitChannels")
        self.tableNodes = self.findChild(QTableWidget, "tableWidget_Nodes")
        self.tableNodes.setColumnWidth(0,40); self.tableNodes.setColumnWidth(1,40)
        self.tableNodes.setColumnWidth(2,40); self.tableNodes.setColumnWidth(3,100)
        self.tableNodes.setColumnWidth(4,100); self.tableNodes.setColumnWidth(5,100)

        # Preprocessing
        self.checkBoxDetrend = self.findChild(QCheckBox, "checkBox_Detrend")
        self.checkBoxDecimation = self.findChild(QCheckBox, "checkBox_Decimation")
        self.lineEditDecimationFactor = self.findChild(QLineEdit, "lineEdit_DecimationFactor")
        self.lineEditInsertPeak = self.findChild(QLineEdit, "lineEdit_InsertPeak")
        self.lineEditSamplingFrequency = self.findChild(QLineEdit, "lineEdit_SamplingFrequency")
        self.listIdentifiedPeaks = self.findChild(QListWidget, "listWidget_DisplayIdentifiedPeaks")
        self.buttonAddIdentifiedPeak = self.findChild(QPushButton, "pushButton_AddIdentifiedPeak")
        self.buttonClearIdentifiedPeaks = self.findChild(QPushButton, "pushButton_ClearIdentifiedPeaks")
        self.buttonDeleteIdentifiedPeak = self.findChild(QPushButton, "pushButton_DeleteIdentifiedPeak")
        self.buttonFDDSvp = self.findChild(QPushButton, "pushButton_FddSvp")
        self.buttonSubmitIdentifiedPeaks = self.findChild(QPushButton, "pushButton_SubmitIdentifiedPeaks")

        # FDD
        self.checkBox_EFDD = self.findChild(QCheckBox, "checkBox_Efdd")
        self.checkBox_FSDD = self.findChild(QCheckBox, "checkBox_Fsdd")
        self.checkBox_OriginalFDD = self.findChild(QCheckBox, "checkBox_OriginalFdd")
        self.comboBoxFDDFigures1 = self.findChild(QComboBox, "comboBox_FddFigures1")
        self.comboBoxFDDFigures2 = self.findChild(QComboBox, "comboBox_FddFigures2")
        self.buttonRunFDD = self.findChild(QPushButton, "pushButton_RunFdd")

        # FDD results
        self.tableDampEfdd = self.findChild(QTableWidget, "tableWidget_DampEfdd")
        self.tableDampFsdd = self.findChild(QTableWidget, "tableWidget_DampFsdd")
        self.tableFreqEfdd = self.findChild(QTableWidget, "tableWidget_FreqEfdd")
        self.tableFreqFdd = self.findChild(QTableWidget, "tableWidget_FreqFdd")
        self.tableFreqFsdd = self.findChild(QTableWidget, "tableWidget_FreqFsdd")
        self.tableModeEfdd = self.findChild(QTableWidget, "tableWidget_ModeEfdd")
        self.tableModeFdd = self.findChild(QTableWidget, "tableWidget_ModeFdd")
        self.tableModeFsdd = self.findChild(QTableWidget, "tableWidget_ModeFsdd")

        # FDD geometry
        self.comboBoxFDDFiguresGeom1 = self.findChild(QComboBox, "comboBox_FddFiguresGeom1")
        self.comboBoxFDDFiguresGeom2 = self.findChild(QComboBox, "comboBox_FddFiguresGeom2")
        self.comboBoxFDDFiguresGeom3 = self.findChild(QComboBox, "comboBox_FddFiguresGeom3")
        self.checkBoxDeformedShapeFDD = self.findChild(QCheckBox, "checkBox_DeformedShapeFDD")
        self.checkBoxValuesOnPlotFDD = self.findChild(QCheckBox, "checkBox_ValuesOnPlotFDD")
        self.lineEditArrowScaleFactorFDD = self.findChild(QLineEdit, "lineEdit_ArrowScaleFactorFDD")
        self.lineEditDimensionScaleFactorXFDD = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_xFDD")
        self.lineEditDimensionScaleFactorYFDD = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_yFDD")
        self.lineEditDimensionScaleFactorZFDD = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_zFDD")
        self.tableDeformedValuesFDD = self.findChild(QTableWidget, "tableWidget_DeformedValuesFDD")

        # SSI
        self.checkBox_SSI_cov = self.findChild(QCheckBox, "checkBox_SsiCov")
        self.checkBox_SSI_dat = self.findChild(QCheckBox, "checkBox_SsiDat")
        self.comboBoxSSIFigures = self.findChild(QComboBox, "comboBox_SsiFigures")
        self.buttonRunSSI = self.findChild(QPushButton, "pushButton_RunSsi")
        # self.lineMinOrder = self.findChild(QLineEdit, "lineEdit_MinOrder")
        self.lineMaxOrder = self.findChild(QLineEdit, "lineEdit_MaxOrder")
        self.lineEditTimeShifts = self.findChild(QLineEdit, "lineEdit_TimeShifts")
        self.lineEditLim0 = self.findChild(QLineEdit, "lineEdit_Lim0")
        self.lineEditLim1 = self.findChild(QLineEdit, "lineEdit_Lim1")
        self.lineEditLim2 = self.findChild(QLineEdit, "lineEdit_Lim2")
        self.lineEditLim3 = self.findChild(QLineEdit, "lineEdit_Lim3")

        # SSI results
        self.tableDampSsiCov = self.findChild(QTableWidget, "tableWidget_DampSsiCov")
        self.tableDampSsiDat = self.findChild(QTableWidget, "tableWidget_DampSsiDat")
        self.tableFreqSsiCov = self.findChild(QTableWidget, "tableWidget_FreqSsiCov")
        self.tableFreqSsiDat = self.findChild(QTableWidget, "tableWidget_FreqSsiDat")
        self.tableModeSsiCov = self.findChild(QTableWidget, "tableWidget_ModeSsiCov")
        self.tableModeSsiDat = self.findChild(QTableWidget, "tableWidget_ModeSsiDat")

        # SSI geometry
        self.comboBoxSSIFiguresGeom1 = self.findChild(QComboBox, "comboBox_SsiFiguresGeom1")
        self.comboBoxSSIFiguresGeom2 = self.findChild(QComboBox, "comboBox_SsiFiguresGeom2")
        self.checkBoxDeformedShapeSSI = self.findChild(QCheckBox, "checkBox_DeformedShapeSSI")
        self.checkBoxValuesOnPlotSSI = self.findChild(QCheckBox, "checkBox_ValuesOnPlotSSI")
        self.lineEditArrowScaleFactorSSI = self.findChild(QLineEdit, "lineEdit_ArrowScaleFactorSSI")
        self.lineEditDimensionScaleFactorXSSI = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_xSSI")
        self.lineEditDimensionScaleFactorYSSI = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_ySSI")
        self.lineEditDimensionScaleFactorZSSI = self.findChild(QLineEdit, "lineEdit_DimensionScaleFactor_zSSI")
        self.tableDeformedValuesSSI = self.findChild(QTableWidget, "tableWidget_DeformedValuesSSI")

        # Action
        # Import Data
        self.buttonCreateFolder.clicked.connect(self.clicker_create_folder)
        self.buttonClearAll.clicked.connect(self.clicker_clear_all)
        self.buttonClearAll.setEnabled(False)
        self.buttonLoadData.clicked.connect(self.clicker_load_data)
        self.buttonLoadData.setEnabled(False)
        self.buttonSubmitSetup.clicked.connect(self.clicker_submit_setup)
        self.buttonSubmitSetup.setEnabled(False)

        # Geometry
        self.buttonLoadNodes.clicked.connect(self.clicker_load_nodes)
        self.buttonLoadNodes.setEnabled(False)
        self.buttonLoadConnectivity.clicked.connect(self.clicker_load_connectivity)
        self.buttonLoadConnectivity.setEnabled(False)
        self.buttonCreateGeometry.clicked.connect(self.clicker_create_geometry)
        self.buttonCreateGeometry.setEnabled(False)
        self.buttonClearChannels.clicked.connect(self.clicker_clear_channels)
        self.buttonSubmitChannels.clicked.connect(self.clicker_submit_channels)
        self.buttonSubmitChannels.setEnabled(False)

        # Pre-proc
        self.buttonAddIdentifiedPeak.clicked.connect(self.add_peak)
        self.buttonDeleteIdentifiedPeak.clicked.connect(self.delete_peak)
        self.buttonClearIdentifiedPeaks.clicked.connect(self.clear_peak)
        self.buttonSubmitIdentifiedPeaks.clicked.connect(self.submit_peak)
        self.buttonAddIdentifiedPeak.setEnabled(False)
        self.buttonDeleteIdentifiedPeak.setEnabled(False)
        self.buttonClearIdentifiedPeaks.setEnabled(False)
        self.buttonSubmitIdentifiedPeaks.setEnabled(False)
        self.buttonFDDSvp.clicked.connect(self.clicker_run_fdd_svp)
        self.buttonFDDSvp.setEnabled(False)

        # FDD
        self.buttonRunFDD.clicked.connect(self.clicker_run_fdd)
        self.buttonRunFDD.setEnabled(False)
        self.comboBoxFDDFigures2.activated.connect(self.plot_fig_fdd)

        # FDD geometry
        self.comboBoxFDDFiguresGeom2.activated.connect(self.plot_fig_fdd_geom)

        # SSI
        self.buttonRunSSI.clicked.connect(self.clicker_run_ssi)
        self.buttonRunSSI.setEnabled(False)
        self.comboBoxSSIFigures.activated.connect(self.plot_fig_ssi)

        # SSI geometry
        self.comboBoxSSIFiguresGeom2.activated.connect(self.plot_fig_ssi_geom)

        self.show()

    # Clicker for selecting the working folder
    def clicker_create_folder(self):
        # Clear previous run data
        self.labelOpenFile.clear()
        while self.tableInputData.rowCount() > 0:
            self.tableInputData.removeRow(0)
        dialog = QFileDialog(self)
        if self.labelCreateFolder.text():
            dialog.setDirectory(self.labelCreateFolder.text())
        dialog.setFileMode(dialog.Directory)

        # we cannot use the native dialog, because we need control over the UI
        options = dialog.Options(dialog.DontUseNativeDialog | dialog.ShowDirsOnly)
        dialog.setOptions(options)

        def check_line_edit(_path):
            if not _path:
                return
            if _path.endswith(QDir.separator()):
                return check_line_edit(_path.rstrip(QDir.separator()))
            _path = QFileInfo(_path)
            if _path.exists() or QFileInfo(_path.absolutePath()).exists():
                button.setEnabled(True)
                return True

        # get the "Open" button in the dialog
        button = dialog.findChild(QDialogButtonBox).button(
            QDialogButtonBox.Open)

        # get the line edit used for the path
        line_edit = dialog.findChild(QLineEdit)
        line_edit.textChanged.connect(check_line_edit)

        # override the existing accept() method, otherwise selectedFiles() will
        # complain about selecting a non-existing path
        def accept():
            if check_line_edit(line_edit.text()):
                # if the path is acceptable, call the base accept() implementation
                QDialog.accept(dialog)
        dialog.accept = accept
        path = None
        if dialog.exec_() and dialog.selectedFiles():
            path = QFileInfo(dialog.selectedFiles()[0]).absoluteFilePath()
            self.labelCreateFolder.setText(path)
        if path:
            self.pathFolder = path
            self.buttonLoadData.setEnabled(True)

    # Clicker for loading the initial data
    def clicker_load_data(self):
        self.tableInputData.setRowCount(0)
        self.tableInputData.setColumnCount(0)
        f_name, f_type = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*);;Text Document (*.txt)"
                                                                            ";;CSV File(*.csv);;Excel Binary File"
                                                                            "(*.xls);;Excel File(*.xlsx)")
        root, extension = os.path.splitext(f_name)
        extent = extension.casefold()
        self.pathLoadFile = root + extent
        if len(os.listdir(self.pathFolder)) != 0:
            self.clicker_clear_all()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Folder Error")
            msg.setText(f'Please select another folder')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        else:
            if f_name:
                shutil.copy(self.pathLoadFile, self.pathFolder)
                self.labelOpenFile.setText(f_name)
                if extent == ".xls" or extent == ".xlsx":
                    data = pd.read_excel(self.pathLoadFile, header=None, sep='\s+', index_col=False, dtype=float)
                elif extent == ".txt" or extent == ".csv":
                    data = pd.read_csv(self.pathLoadFile, header=None, sep='\s+', index_col=False, dtype=float)
                else:
                    pass
                data = data.to_numpy()
                time_history = data.shape[0]
                no_channels = data.shape[1]
                self.tableInputData.setColumnCount(no_channels)
                self.tableInputData.setRowCount(6)
                for i in list(np.arange(0,6)):
                    for j in range(no_channels):
                        if i == 5:
                            self.tableInputData.setItem(i, j, QTableWidgetItem('...'))
                        else:
                            self.tableInputData.setItem(i, j, QTableWidgetItem(str(data[i,j])))
                self.labelTimeSteps.setText(str(data.shape[0]))
                self.labelNumOfChannels.setText(str(data.shape[1]))
                self.inputData = data
                self.no_channels = no_channels
                headers_labels = []
                for i in range(no_channels):
                    _temp = f'Channel_{i + 1}'
                    headers_labels.append(_temp)
                self.tableInputData.setHorizontalHeaderLabels(headers_labels)
                self.tableInputData.horizontalHeader().sectionDoubleClicked.connect(self.change_horizontal_header)
                self.num_setup = '1'
                self.channelNamesDict[self.num_setup] = headers_labels
            QMessageBox.about(self, "Selected File", "Successfully loaded!")
            self.buttonSubmitSetup.setEnabled(True)
            self.buttonClearAll.setEnabled(True)

    # Change headers of load data
    def change_horizontal_header(self, index):
        old_header = self.tableInputData.horizontalHeaderItem(index).text()
        new_header, ok = QInputDialog.getText(self,
                                                   'Change header label for column %d' % index,
                                                   'Header:',
                                                   QLineEdit.Normal,
                                                   old_header)
        if ok:
            self.tableInputData.horizontalHeaderItem(index).setText(new_header)

    # Submit number of Set up
    def clicker_submit_setup(self):
        headers_labels = []
        for c in range(self.tableInputData.columnCount()):
            it = self.tableInputData.horizontalHeaderItem(c)
            dict_headers = it.text()
            headers_labels.append(dict_headers)
        shutil.copy(self.pathLoadFile, self.pathFolder)
        results_directory = self.pathFolder + '/Results'
        os.makedirs(results_directory)
        self.resultsDirectory = results_directory
        self.channelNamesDict[self.num_setup] = headers_labels
        self.buttonFDDSvp.setEnabled(True)
        self.buttonSubmitChannels.setEnabled(True)
        self.comboBoxChannelNames.clear()
        for item in headers_labels:
            self.comboBoxChannelNames.addItem(item)
        self.buttonSubmitSetup.setEnabled(False)
        self.buttonLoadConnectivity.setEnabled(True)
        self.buttonLoadNodes.setEnabled(True)

    # Clicker for loading the nodes' coordinates
    def clicker_load_nodes(self):
        f_name, f_type = QFileDialog.getOpenFileName(self, "Open File", "",
                                                     "All Files (*);;Text Document (*.txt)"
                                                     ";;CSV File(*.csv);;Excel Binary File"
                                                     "(*.xls);;Excel File(*.xlsx)")
        root, extension = os.path.splitext(f_name)
        extent = extension.casefold()
        path_load_nodes = root + extent
        if f_name:
            self.labelOpenFileNodes.setText(f_name)
            if extent in [".xls", ".xlsx"]:
                nodes = pd.read_excel(path_load_nodes, header=None, sep='\s+', index_col=False, dtype=float)
            elif extent in [".txt", ".csv"]:
                nodes = pd.read_csv(path_load_nodes, header=None, sep='\s+', index_col=False, dtype=float)
            else:
                pass
            self.nodes = nodes.to_numpy()
            num_nodes = self.nodes.shape[0]
            self.tableNodes.setRowCount(num_nodes)
            for i in range(num_nodes):
                p = Point()
                p.xyz = [self.nodes[i,0], self.nodes[i,1], self.nodes[i,2]]
                p.xyz_new = [self.nodes[i,0], self.nodes[i,1], self.nodes[i,2]]
                self.nodesDict[i + 1] = p
                self.tableNodes.setItem(i, 0, QTableWidgetItem(str(p.xyz[0])))
                self.tableNodes.setItem(i, 1, QTableWidgetItem(str(p.xyz[1])))
                self.tableNodes.setItem(i, 2, QTableWidgetItem(str(p.xyz[2])))
            self.changed_items = []
            self.changed_items_name = []
            self.tableNodes.itemChanged.connect(self.log_change)
            self.setup_clicks = 2
            if len(self.connectivity) !=0:
                self.buttonCreateGeometry.setEnabled(True)
            self.buttonSubmitChannels.setEnabled(True)

    # Trace the changes in table
    def log_change(self, item):
        self.changed_items.append(item)
        self.changed_items_name.append(item.text())

    # Clicker for loading the connectivity of the structure
    def clicker_load_connectivity(self):
        f_name, f_type = QFileDialog.getOpenFileName(self, "Open File", "",
                                                     "All Files (*);;Text Document (*.txt)"
                                                     ";;CSV File(*.csv);;Excel Binary File"
                                                     "(*.xls);;Excel File(*.xlsx)")
        root, extension = os.path.splitext(f_name)
        extent = extension.casefold()
        path_load_connectivity = root + extent
        if f_name:
            self.labelOpenFileConnectivity.setText(f_name)
            if extent in [".xls", ".xlsx"]:
                connectivity = pd.read_excel(path_load_connectivity, header=None, sep='\s+', index_col=False, dtype=int)
            elif extent in [".txt", ".csv"]:
                connectivity = pd.read_csv(path_load_connectivity, header=None, sep='\s+', index_col=False, dtype=int)
            else:
                pass
            self.connectivity = connectivity.to_numpy()
            if len(self.nodes) !=0:
                self.buttonCreateGeometry.setEnabled(True)

    # Clicker for creating the initial geometry
    def clicker_create_geometry(self):
        connectivity = self.connectivity
        nodes = self.nodes
        num_frames = connectivity.shape[0]
        setattr(Axes3D, 'annotate3D', drawing_tools_3d.annotate3d)
        _fig = plt.figure()
        ax = plt.axes(projection="3d")
        for k in range(num_frames):
            x1 = nodes[connectivity[k, 0] - 1, 0]
            y1 = nodes[connectivity[k, 0] - 1, 1]
            z1 = nodes[connectivity[k, 0] - 1, 2]
            x2 = nodes[connectivity[k, 1] - 1, 0]
            y2 = nodes[connectivity[k, 1] - 1, 1]
            z2 = nodes[connectivity[k, 1] - 1, 2]
            xx = [x1, x2]; yy = [y1, y2]; zz = [z1, z2]
            ax.plot3D(xx, yy, zz, 'black')
        for i in range(nodes.shape[0]):
            xs = nodes[i, 0]; ys = nodes[i, 1]; zs = nodes[i, 2]
            ax.scatter(xs, ys, zs, color="black", marker='o')
            ax.annotate3D(f'P{i + 1}', (xs, ys, zs), xytext=(3, 3), textcoords='offset points')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        x = ax.get_xlim3d()
        y = ax.get_ylim3d()
        z = ax.get_zlim3d()
        a = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]
        b = np.amax(a)
        ax.set_xlim3d(x[0] - (b - a[0]) / 2, x[1] + (b - a[0]) / 2)
        ax.set_ylim3d(y[0] - (b - a[1]) / 2, y[1] + (b - a[1]) / 2)
        ax.set_zlim3d(z[0] - (b - a[2]) / 2, z[1] + (b - a[2]) / 2)
        # Add layouts + widgets
        self.plotWidgetGeometry = FigureCanvas(_fig)
        self.layGeometryPlot.addWidget(self.plotWidgetGeometry)
        _toolWidget = NavigationToolbar(self.plotWidgetGeometry, self)
        self.layGeometryTool.addWidget(_toolWidget)
        self.buttonCreateGeometry.setEnabled(False)
        plt.close('all')

    # Clicker for submit channels in DOF
    def clicker_submit_channels(self):
        _channelNames = self.channelNamesDict[self.num_setup]
        check = all(elem in self.changed_items_name for elem in _channelNames)
        if not check:
            if len(self.nodes) != 0:
                for key, value in self.nodesDict.items():
                    self.nodesDict[key].channels = []
                self.tableNodes.itemChanged.disconnect(self.log_change)
                while self.tableNodes.rowCount() > 0:
                    self.tableNodes.removeRow(0)
                num_nodes = self.nodes.shape[0]
                self.tableNodes.setRowCount(num_nodes)
                for i in range(num_nodes):
                    p = self.nodesDict[i + 1]
                    self.tableNodes.setItem(i, 0, QTableWidgetItem(str(p.xyz[0])))
                    self.tableNodes.setItem(i, 1, QTableWidgetItem(str(p.xyz[1])))
                    self.tableNodes.setItem(i, 2, QTableWidgetItem(str(p.xyz[2])))
                self.changed_items = []
                self.changed_items_name = []
                self.tableNodes.itemChanged.connect(self.log_change)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Channels Name Error")
            msg.setText(f'Please assign the channels to the equivalent DOF')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        elif len(self.changed_items) != self.no_channels:
            if len(self.nodes) != 0:
                for key, value in self.nodesDict.items():
                    self.nodesDict[key].channels = []
                self.tableNodes.itemChanged.disconnect(self.log_change)
                while self.tableNodes.rowCount() > 0:
                    self.tableNodes.removeRow(0)
                num_nodes = self.nodes.shape[0]
                self.tableNodes.setRowCount(num_nodes)
                for i in range(num_nodes):
                    p = self.nodesDict[i + 1]
                    self.tableNodes.setItem(i, 0, QTableWidgetItem(str(p.xyz[0])))
                    self.tableNodes.setItem(i, 1, QTableWidgetItem(str(p.xyz[1])))
                    self.tableNodes.setItem(i, 2, QTableWidgetItem(str(p.xyz[2])))
                self.changed_items = []
                self.changed_items_name = []
                self.tableNodes.itemChanged.connect(self.log_change)
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Channels Number Error")
            msg.setText(f'Please assign {self.no_channels} channels to the equivalent DOF')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        else:
            for i in range(len(self.changed_items)):
                _node = self.changed_items[i]
                p = self.nodesDict[_node.row() + 1]
                if len(p.channels) == 0:
                    p.channels = [0, 0, 0]
                p.channels[_node.column() - 3] = _node.text()
                self.channelsDict[_node.text()] = [_node.row() + 1, _node.column() - 3]
            self.buttonSubmitChannels.setEnabled(False)
            if len(self.channelsDict) != 0:
                self.comboBoxChannelAssignments.clear()
                _channelNames = self.channelNamesDict[self.num_setup]
                for _name in _channelNames:
                    temp = self.channelsDict[_name]
                    id_node = temp[0]
                    if temp[1] == 0:
                        temp_string = f'The {_name} is assigned in Node {id_node} in x direction'
                        self.comboBoxChannelAssignments.addItem(temp_string)
                    elif temp[1] == 1:
                        temp_string = f'The {_name} is assigned in Node {id_node} in y direction'
                        self.comboBoxChannelAssignments.addItem(temp_string)
                    elif temp[1] == 2:
                        temp_string = f'The {_name} is assigned in Node {id_node} in z direction'
                        self.comboBoxChannelAssignments.addItem(temp_string)

    # Clicker to clear channels' assignments in DOF
    def clicker_clear_channels(self):
        if len(self.nodes) !=0:
            for key, value in self.nodesDict.items():
                self.nodesDict[key].channels = []
            self.tableNodes.itemChanged.disconnect(self.log_change)
            while self.tableNodes.rowCount() > 0:
                self.tableNodes.removeRow(0)
            num_nodes = self.nodes.shape[0]
            self.tableNodes.setRowCount(num_nodes)
            for i in range(num_nodes):
                p = self.nodesDict[i + 1]
                self.tableNodes.setItem(i, 0, QTableWidgetItem(str(p.xyz[0])))
                self.tableNodes.setItem(i, 1, QTableWidgetItem(str(p.xyz[1])))
                self.tableNodes.setItem(i, 2, QTableWidgetItem(str(p.xyz[2])))
            self.changed_items = []
            self.changed_items_name = []
            self.tableNodes.itemChanged.connect(self.log_change)
            self.buttonSubmitChannels.setEnabled(True)

    # Clicker for running the FDD_SVP
    def clicker_run_fdd_svp(self):
        data = self.inputData
        self.parameters = Parameters()
        self.parameters.sampling_frequency = float(self.lineEditSamplingFrequency.text())
        self.parameters.decimation_factor = int(self.lineEditDecimationFactor.text())
        if self.checkBoxDetrend.isChecked():
            data = signal.detrend(data, axis=0) # Trend removal
        if self.checkBoxDecimation.isChecked():
            q = self.parameters.decimation_factor  # Decimation factor
            data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
            self.parameters.sampling_frequency = self.parameters.sampling_frequency/q  # [Hz] Decimated sampling freq.
        self.inputData = data
        # Run FDD
        fdd = oma.FDDsvp(self.inputData,  self.parameters.sampling_frequency)
        fdd[0].savefig(self.resultsDirectory + '/' + 'SV(PSD)_plot.png')
        self.plotWidgetFddSvp = FigureCanvas(fdd[0])
        self.layFddSvpPlot.addWidget(self.plotWidgetFddSvp)
        _toolWidget = NavigationToolbar(self.plotWidgetFddSvp, self)
        self.layFddSvpTool.addWidget(_toolWidget)
        self.fddhelper = fdd[1]
        self.buttonAddIdentifiedPeak.setEnabled(True)
        self.buttonDeleteIdentifiedPeak.setEnabled(True)
        self.buttonClearIdentifiedPeaks.setEnabled(True)
        self.buttonSubmitIdentifiedPeaks.setEnabled(True)
        self.buttonFDDSvp.setEnabled(False)
        cursor = mplcursors.cursor(hover=False)
        @cursor.connect("add")
        def cursor_clicked(sel):
            global my_clicker_helper
            my_clicker_helper = sel.target[0]

    # Clicker for adding an identified peak
    def add_peak(self):
        item = self.lineEditInsertPeak.text()
        if item == "":
             global my_clicker_helper
             self.listIdentifiedPeaks.addItem("{:.2f}".format(my_clicker_helper))
             "{:.6f}".format(my_clicker_helper)
             my_clicker_helper = 0
        else:
             self.listIdentifiedPeaks.addItem(item)
             self.lineEditInsertPeak.setText("")

    # Clicker for deleting an identified peak
    def delete_peak(self):
        clicked = self.listIdentifiedPeaks.currentRow()
        self.listIdentifiedPeaks.takeItem(clicked)

    # Clicker for clear the list
    def clear_peak(self):
        self.listIdentifiedPeaks.clear()
        self.buttonSubmitIdentifiedPeaks.setEnabled(True)

    # Clicker for submit the identified peaks
    def submit_peak(self):
        num = self.listIdentifiedPeaks.count()
        if num == 0:
            self.clear_peak()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Identified Peaks Number Error")
            msg.setText(f'Please assign at least one identified peaks to continue')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
        else:
            self.parameters.identified_peaks = []
            for i in range(num):
                self.parameters.identified_peaks.append(float(self.listIdentifiedPeaks.item(i).text()))
            self.buttonSubmitIdentifiedPeaks.setEnabled(False)
            self.buttonRunFDD.setEnabled(True)
            self.buttonRunSSI.setEnabled(True)
            plt.close('all')

    # Clicker for running the FDD
    def clicker_run_fdd(self):
        results_directory = self.resultsDirectory + '/FDD'
        os.makedirs(results_directory)
        self.comboBoxFDDFiguresGeom3.addItem(self.num_setup)
        for i in range(self.listIdentifiedPeaks.count()):
            item = f'Mode{i + 1}'
            self.comboBoxFDDFigures2.addItem(item)
            self.comboBoxFDDFiguresGeom2.addItem(item)
        if self.checkBox_OriginalFDD.isChecked():
            res_fdd = oma.FDDmodEX(self.parameters.identified_peaks, self.fddhelper)
            self.write_to_txt(res_fdd, results_directory, 'FDD')
            self.write_to_gui(res_fdd, self.tableFreqFdd, 'Frequencies', 'FDD')
            self.write_to_gui(res_fdd, self.tableModeFdd, 'Mode Shapes', 'FDD')
            self.comboBoxFDDFiguresGeom1.addItem('FDD')
            for i in range(self.listIdentifiedPeaks.count()):
                _modes = res_fdd['Mode Shapes'].real
                modes = _modes.T
                _temp = f'Mode{i + 1}'
                self.modesDict[self.num_setup + 'FDD' + _temp] = modes[i]
        if self.checkBox_EFDD.isChecked():
            _fig, res_efdd = oma.EFDDmodEX(self.parameters.identified_peaks, self.fddhelper, method='EFDD', plot=True)
            self.write_to_txt(res_efdd, results_directory, 'EFDD')
            self.write_to_gui(res_efdd, self.tableFreqEfdd, 'Frequencies', 'EFDD')
            self.write_to_gui(res_efdd, self.tableDampEfdd, 'Damping', 'EFDD')
            self.write_to_gui(res_efdd, self.tableModeEfdd, 'Mode Shapes', 'EFDD')
            self.comboBoxFDDFigures1.addItem('EFDD')
            self.comboBoxFDDFiguresGeom1.addItem('EFDD')
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'Mode{i + 1}'
                self._figuresFDD['EFDD' + _temp] = _fig[i]
                _modes = res_efdd['Mode Shapes'].real
                modes = _modes.T
                self.modesDict[self.num_setup + 'EFDD' + _temp] = modes[i]
                _temp_png = f'EFDDMode{i + 1}.png'
                _fig[i].savefig(results_directory + '/' + _temp_png)
            self.plotWidgetFdd = FigureCanvas(_fig[0])
        if self.checkBox_FSDD.isChecked():
            _fig, res_fsdd = oma.EFDDmodEX(self.parameters.identified_peaks, self.fddhelper, method='FSDD', npmax = 35,
                                           MAClim=0.95, plot=True)
            self.write_to_txt(res_fsdd, results_directory, 'FSDD')
            self.write_to_gui(res_fsdd, self.tableFreqFsdd, 'Frequencies', 'FSDD')
            self.write_to_gui(res_fsdd, self.tableDampFsdd, 'Damping', 'FSDD')
            self.write_to_gui(res_fsdd, self.tableModeFsdd, 'Mode Shapes', 'FSDD')
            self.comboBoxFDDFigures1.addItem('FSDD')
            self.comboBoxFDDFiguresGeom1.addItem('FSDD')
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'Mode{i + 1}'
                self._figuresFDD['FSDD' + _temp] = _fig[i]
                _modes = res_fsdd['Mode Shapes'].real
                modes = _modes.T
                self.modesDict[self.num_setup + 'FSDD' + _temp] = modes[i]
                _temp_png = f'FSDDMode{i + 1}.png'
                _fig[i].savefig(results_directory + '/' + _temp_png)
            self.plotWidgetFdd = FigureCanvas(_fig[0])
        self.fddhelper = 0
        # Add layouts + widgets
        # self.plotWidgetFdd = FigureCanvas()
        self.layFddPlot.addWidget(self.plotWidgetFdd)
        _toolWidget = NavigationToolbar(self.plotWidgetFdd, self)
        self.layFddTool.addWidget(_toolWidget)
        self.buttonRunFDD.setEnabled(False)
        plt.close('all')

    # Clicker to display FDD results
    def plot_fig_fdd(self):
        self.remove_widget(self.layFddTool)
        method_name = self.comboBoxFDDFigures1.currentText()
        mode_name = self.comboBoxFDDFigures2.currentText()
        self.layFddPlot.removeWidget(self.plotWidgetFdd)
        _fig = self._figuresFDD.get(method_name + mode_name)
        self.plotWidgetFdd = FigureCanvas(_fig)
        self.layFddPlot.addWidget(self.plotWidgetFdd)
        self.layFddTool.addWidget(NavigationToolbar(self.plotWidgetFdd, self))
        plt.close(_fig)

    # Clicker to display FDD geometry
    def plot_fig_fdd_geom(self): #modificata la function
        if self.lineEditDimensionScaleFactorXFDD.text().isdigit() and int(self.lineEditDimensionScaleFactorXFDD.text()) > 0\
        and self.lineEditDimensionScaleFactorYFDD.text().isdigit() and int(self.lineEditDimensionScaleFactorYFDD.text()) > 0\
        and self.lineEditDimensionScaleFactorZFDD.text().isdigit() and int(self.lineEditDimensionScaleFactorZFDD.text()) > 0\
        and self.lineEditArrowScaleFactorFDD.text().isdigit() and int(self.lineEditArrowScaleFactorFDD.text()) > 0:
            setattr(Axes3D, 'annotate3D', drawing_tools_3d.annotate3d)
            setattr(Axes3D, 'arrow3D', drawing_tools_3d.arrow3d)
            self.remove_widget(self.layFddToolGeom)
            self.layFddPlotGeom.removeWidget(self.plotWidgetFddGeom)
            method_name = self.comboBoxFDDFiguresGeom1.currentText()
            mode_name = self.comboBoxFDDFiguresGeom2.currentText()
            num_setup = self.comboBoxFDDFiguresGeom3.currentText()
            _mode = self.modesDict[num_setup + method_name + mode_name]
            _sum = 0
            _channelNames = self.channelNamesDict[num_setup]
            scale_factor_length_x = int(self.lineEditDimensionScaleFactorXFDD.text())
            scale_factor_length_y = int(self.lineEditDimensionScaleFactorYFDD.text())
            scale_factor_length_z = int(self.lineEditDimensionScaleFactorZFDD.text())
            for channel in _channelNames:
                factor = 1
                if channel[0] == '-':
                    factor = -1
                node = self.channelsDict[channel]
                point = self.nodesDict[node[0]]
                if node[1] == 0:
                    point.xyz_new[0] += scale_factor_length_x * _mode[_sum]
                elif node[1] == 1:
                    point.xyz_new[1] += scale_factor_length_y * _mode[_sum]
                elif node[1] == 2:
                    point.xyz_new[2] += scale_factor_length_z * _mode[_sum]
                _sum += 1
            connectivity = self.connectivity
            nodes = self.nodes
            num_frames = connectivity.shape[0]
            _fig = plt.figure()
            ax = plt.axes(projection="3d")
            scale_factor_arrow = int(self.lineEditArrowScaleFactorFDD.text())
            for k in range(num_frames):
                p1 = self.nodesDict[connectivity[k, 0]]
                p2 = self.nodesDict[connectivity[k, 1]]
                x1 = p1.xyz[0]
                y1 = p1.xyz[1]
                z1 = p1.xyz[2]
                x2 = p2.xyz[0]
                y2 = p2.xyz[1]
                z2 = p2.xyz[2]
                xx = [x1, x2]
                yy = [y1, y2]
                zz = [z1, z2]
                ax.plot3D(xx, yy, zz, 'black')
                if self.checkBoxDeformedShapeFDD.isChecked():
                    x1 = p1.xyz_new[0]
                    y1 = p1.xyz_new[1]
                    z1 = p1.xyz_new[2]
                    x2 = p2.xyz_new[0]
                    y2 = p2.xyz_new[1]
                    z2 = p2.xyz_new[2]
                    xx = [x1, x2]
                    yy = [y1, y2]
                    zz = [z1, z2]
                    ax.plot3D(xx, yy, zz, 'red')
            self.tableDeformedValuesFDD.setRowCount(nodes.shape[0])
            for i in range(nodes.shape[0]):
                p1 = self.nodesDict[i + 1]
                xs1 = p1.xyz[0]
                ys1 = p1.xyz[1]
                zs1 = p1.xyz[2]
                ax.scatter(xs1, ys1, zs1, color="black", marker='o')
                ax.annotate3D(f'P{i + 1}', (xs1, ys1, zs1), xytext=(3, 3),
                              textcoords='offset points')
                xs2 = p1.xyz_new[0]
                ys2 = p1.xyz_new[1]
                zs2 = p1.xyz_new[2]
                if self.checkBoxDeformedShapeFDD.isChecked():
                    ax.scatter(xs2, ys2, zs2, color="red", marker='o')
                _x = round((xs2 - xs1), 2)
                _y = round((ys2 - ys1), 2)
                _z = round((zs2 - zs1), 2)
                self.tableDeformedValuesFDD.setItem(i, 0, QTableWidgetItem(f'P{i + 1}'))
                self.tableDeformedValuesFDD.setItem(i, 1, QTableWidgetItem(str(_x)))
                self.tableDeformedValuesFDD.setItem(i, 2, QTableWidgetItem(str(_y)))
                self.tableDeformedValuesFDD.setItem(i, 3, QTableWidgetItem(str(_z)))
                if self.checkBoxValuesOnPlotFDD.isChecked():
                    ax.annotate3D(f'P{i+1}: dx={_x}, dy={_y}, dz={_z}',
                                  (xs1, ys1, zs1), xytext=(3, 3), textcoords='offset points')
                else:
                    ax.annotate3D(f'P{i + 1}', (xs1, ys1, zs1), xytext=(3, 3),
                                  textcoords='offset points')
                ax.arrow3D(xs1, ys1, zs1,
                           xs2 - xs1, ys2 - ys1, zs2 - zs1,
                           mutation_scale=7,
                           arrowstyle="-|>",
                           color="r",
                           lw=scale_factor_arrow * 1)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            x = ax.get_xlim3d()
            y = ax.get_ylim3d()
            z = ax.get_zlim3d()
            a = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]
            b = np.amax(a)
            ax.set_xlim3d(x[0] - (b - a[0]) / 2, x[1] + (b - a[0]) / 2)
            ax.set_ylim3d(y[0] - (b - a[1]) / 2, y[1] + (b - a[1]) / 2)
            ax.set_zlim3d(z[0] - (b - a[2]) / 2, z[1] + (b - a[2]) / 2)
            _sum = 0
            _channelNames = self.channelNamesDict[num_setup]
            for channel in _channelNames:
                node = self.channelsDict[channel]
                point = self.nodesDict[node[0]]
                if node[1] == 0:
                    _asd = point.xyz[0]
                    point.xyz_new[0] = _asd
                elif node[1] == 1:
                    _asd = point.xyz[1]
                    point.xyz_new[1] = _asd
                elif node[1] == 2:
                    _asd = point.xyz[2]
                    point.xyz_new[2] = _asd
                _sum += 1
            self.plotWidgetFddGeom = FigureCanvas(_fig)
            self.layFddPlotGeom.addWidget(self.plotWidgetFddGeom)
            self.layFddToolGeom.addWidget(NavigationToolbar(self.plotWidgetFddGeom, self))
            plt.close(_fig)
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Scale Error")
            msg.setText(f'Please enter a positive integer')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

    # Clicker for running the SSI
    def clicker_run_ssi(self):
        self.parameters.time_shifts = int(self.lineEditTimeShifts.text())
        # self.parameters.min_order = int(self.lineMinOrder.text())
        lim0 = float(self.lineEditLim0.text())
        lim1 = float(self.lineEditLim1.text())
        lim2 = float(self.lineEditLim2.text())
        lim3 = float(self.lineEditLim3.text())
        limVal = (lim0, lim1, lim2, lim3)
        if self.lineMaxOrder.text() != '':
            self.parameters.max_order = int(self.lineMaxOrder.text())
        else:
            self.parameters.max_order = None



        # Extract the modal properties
        results_directory = self.resultsDirectory + '/SSI'
        os.makedirs(results_directory)
        for i in range(self.listIdentifiedPeaks.count()):
            item = f'Mode{i + 1}'
            self.comboBoxSSIFiguresGeom2.addItem(item)
        if self.checkBox_SSI_cov.isChecked():
            ssi_cov = oma.SSIcovStaDiag(self.inputData, self.parameters.sampling_frequency, self.parameters.time_shifts,
                                        # ordmin=self.parameters.min_order, 
                                        ordmax=self.parameters.max_order,
                                        lim=limVal)
            self.plotWidgetSsi = FigureCanvas(ssi_cov[0])
            _temp = f'SSI_cov_Results'
            self._figuresSSI[_temp] = ssi_cov[0]
            _temp_png = f'SSI_cov_Results.png'
            ssi_cov[0].savefig(results_directory + '/' + _temp_png)
            res_ssi_cov = oma.SSIModEX(self.parameters.identified_peaks, ssi_cov[1])
            self.write_to_txt(res_ssi_cov, results_directory, 'SSIcov')
            self.write_to_gui(res_ssi_cov, self.tableFreqSsiCov, 'Frequencies', 'SSIcov')
            self.write_to_gui(res_ssi_cov, self.tableDampSsiCov, 'Damping', 'SSIcov')
            self.write_to_gui(res_ssi_cov, self.tableModeSsiCov, 'Mode Shapes', 'SSIcov')
            self.comboBoxSSIFiguresGeom1.addItem('SSIcov')
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'SSIcovMode{i + 1}'
                _modes = res_ssi_cov['Mode Shapes'].real
                modes = _modes.T
                self.modesDict[_temp] = modes[i]
        if self.checkBox_SSI_dat.isChecked():
            ssi_dat = oma.SSIdatStaDiag(self.inputData, self.parameters.sampling_frequency, self.parameters.time_shifts,
                                        ordmax=self.parameters.max_order,
                                        lim=limVal)
            self.plotWidgetSsi = FigureCanvas(ssi_dat[0])
            _temp = f'SSI_dat_Results'
            self._figuresSSI[_temp] = ssi_dat[0]
            _temp_png = f'SSI_dat_Results.png'
            ssi_dat[0].savefig(results_directory + '/' + _temp_png)
            res_ssi_dat = oma.SSIModEX(self.parameters.identified_peaks, ssi_dat[1])
            self.write_to_txt(res_ssi_dat, results_directory, 'SSIdat')
            self.write_to_gui(res_ssi_dat, self.tableFreqSsiDat, 'Frequencies', 'SSIdat')
            self.write_to_gui(res_ssi_dat, self.tableDampSsiDat, 'Damping', 'SSIdat')
            self.write_to_gui(res_ssi_dat, self.tableModeSsiDat, 'Mode Shapes', 'SSIdat')
            self.comboBoxSSIFiguresGeom1.addItem('SSIdat')
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'SSIdatMode{i + 1}'
                _modes = res_ssi_dat['Mode Shapes'].real
                modes = _modes.T
                self.modesDict[_temp] = modes[i]
        _tempList = self._figuresSSI.keys()
        for item in _tempList:
            self.comboBoxSSIFigures.addItem(item)
        # Add layouts + widgets
        # self.plotWidgetSsi = FigureCanvas()
        self.laySsiPlot.addWidget(self.plotWidgetSsi)
        _toolWidget = NavigationToolbar(self.plotWidgetSsi, self)
        self.laySsiTool.addWidget(_toolWidget)
        self.buttonRunSSI.setEnabled(False)
        plt.close('all')

    # Clicker to display SSI results
    def plot_fig_ssi(self):
        self.remove_widget(self.laySsiTool)
        result_name = self.comboBoxSSIFigures.currentText()
        self.laySsiPlot.removeWidget(self.plotWidgetSsi)
        _fig = self._figuresSSI.get(result_name)
        self.plotWidgetSsi = FigureCanvas(_fig)
        self.laySsiPlot.addWidget(self.plotWidgetSsi)
        self.laySsiTool.addWidget(NavigationToolbar(self.plotWidgetSsi, self))
        plt.close(_fig)

    # Clicker to display SSI geometry
    def plot_fig_ssi_geom(self): #modificata la function
        if self.lineEditDimensionScaleFactorXSSI.text().isdigit() and int(self.lineEditDimensionScaleFactorXSSI.text()) > 0\
        and self.lineEditDimensionScaleFactorYSSI.text().isdigit() and int(self.lineEditDimensionScaleFactorYSSI.text()) > 0\
        and self.lineEditDimensionScaleFactorZSSI.text().isdigit() and int(self.lineEditDimensionScaleFactorZSSI.text()) > 0\
        and self.lineEditArrowScaleFactorSSI.text().isdigit() and int(self.lineEditArrowScaleFactorSSI.text()) > 0:
            setattr(Axes3D, 'annotate3D', drawing_tools_3d.annotate3d)
            setattr(Axes3D, 'arrow3D', drawing_tools_3d.arrow3d)
            self.remove_widget(self.laySsiToolGeom)
            self.laySsiPlotGeom.removeWidget(self.plotWidgetSsiGeom)
            method_name = self.comboBoxSSIFiguresGeom1.currentText()
            mode_name = self.comboBoxSSIFiguresGeom2.currentText()
            _mode = self.modesDict[method_name + mode_name]
            _sum = 0
            _channelNames = self.channelNamesDict[self.num_setup]
            scale_factor_length_x = int(self.lineEditDimensionScaleFactorXSSI.text())
            scale_factor_length_y = int(self.lineEditDimensionScaleFactorYSSI.text())
            scale_factor_length_z = int(self.lineEditDimensionScaleFactorZSSI.text())
            for channel in _channelNames:
                factor = 1
                if channel[0] == '-':
                    factor = -1
                node = self.channelsDict[channel]
                point = self.nodesDict[node[0]]
                if node[1] == 0:
                    point.xyz_new[0] += scale_factor_length_x * _mode[_sum]
                elif node[1] == 1:
                    point.xyz_new[1] += scale_factor_length_y * _mode[_sum]
                elif node[1] == 2:
                    point.xyz_new[2] += scale_factor_length_z * _mode[_sum]
                _sum += 1
            connectivity = self.connectivity
            nodes = self.nodes
            num_frames = connectivity.shape[0]
            _fig = plt.figure()
            ax = plt.axes(projection="3d")
            scale_factor_arrow = int(self.lineEditArrowScaleFactorSSI.text())
            for k in range(num_frames):
                p1 = self.nodesDict[connectivity[k, 0]]
                p2 = self.nodesDict[connectivity[k, 1]]
                x1 = p1.xyz[0]
                y1 = p1.xyz[1]
                z1 = p1.xyz[2]
                x2 = p2.xyz[0]
                y2 = p2.xyz[1]
                z2 = p2.xyz[2]
                xx = [x1, x2]
                yy = [y1, y2]
                zz = [z1, z2]
                ax.plot3D(xx, yy, zz, 'black')
                if self.checkBoxDeformedShapeSSI.isChecked():
                    x1 = p1.xyz_new[0]
                    y1 = p1.xyz_new[1]
                    z1 = p1.xyz_new[2]
                    x2 = p2.xyz_new[0]
                    y2 = p2.xyz_new[1]
                    z2 = p2.xyz_new[2]
                    xx = [x1, x2]
                    yy = [y1, y2]
                    zz = [z1, z2]
                    ax.plot3D(xx, yy, zz, 'red')
            self.tableDeformedValuesSSI.setRowCount(nodes.shape[0])
            for i in range(nodes.shape[0]):
                p1 = self.nodesDict[i + 1]
                xs1 = p1.xyz[0]
                ys1 = p1.xyz[1]
                zs1 = p1.xyz[2]
                ax.scatter(xs1, ys1, zs1, color="black", marker='o')
                ax.annotate3D(f'P{i + 1}', (xs1, ys1, zs1), xytext=(3, 3),
                              textcoords='offset points')
                xs2 =  p1.xyz_new[0]
                ys2 =  p1.xyz_new[1]
                zs2 =  p1.xyz_new[2]
                if self.checkBoxDeformedShapeSSI.isChecked():
                    ax.scatter(xs2, ys2, zs2, color="red", marker='o')
                _x = round((xs2 - xs1), 2)
                _y = round((ys2 - ys1), 2)
                _z = round((zs2 - zs1), 2)
                self.tableDeformedValuesSSI.setItem(i, 0, QTableWidgetItem(f'P{i + 1}'))
                self.tableDeformedValuesSSI.setItem(i, 1, QTableWidgetItem(str(_x)))
                self.tableDeformedValuesSSI.setItem(i, 2, QTableWidgetItem(str(_y)))
                self.tableDeformedValuesSSI.setItem(i, 3, QTableWidgetItem(str(_z)))
                if self.checkBoxValuesOnPlotSSI.isChecked():
                    ax.annotate3D(f'P{i+1}: dx={_x}, dy={_y}, dz={_z}',
                                  (xs1, ys1, zs1), xytext=(3, 3), textcoords='offset points')
                else:
                    ax.annotate3D(f'P{i + 1}', (xs1, ys1, zs1), xytext=(3, 3),
                                  textcoords='offset points')
                ax.arrow3D(xs1, ys1, zs1,
                           xs2 - xs1, ys2 - ys1, zs2 - zs1,
                           mutation_scale=7,
                           arrowstyle="-|>",
                           color="r",
                           lw=scale_factor_arrow * 1)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            x = ax.get_xlim3d()
            y = ax.get_ylim3d()
            z = ax.get_zlim3d()
            a = [x[1] - x[0], y[1] - y[0], z[1] - z[0]]
            b = np.amax(a)
            ax.set_xlim3d(x[0] - (b - a[0]) / 2, x[1] + (b - a[0]) / 2)
            ax.set_ylim3d(y[0] - (b - a[1]) / 2, y[1] + (b - a[1]) / 2)
            ax.set_zlim3d(z[0] - (b - a[2]) / 2, z[1] + (b - a[2]) / 2)
            _sum = 0
            _channelNames = self.channelNamesDict[self.num_setup]
            for channel in _channelNames:
                node = self.channelsDict[channel]
                point = self.nodesDict[node[0]]
                if node[1] == 0:
                    _asd = point.xyz[0]
                    point.xyz_new[0] = _asd
                elif node[1] == 1:
                    _asd = point.xyz[1]
                    point.xyz_new[1] = _asd
                elif node[1] == 2:
                    _asd = point.xyz[2]
                    point.xyz_new[2] = _asd
                _sum += 1
            self.plotWidgetSsiGeom = FigureCanvas(_fig)
            self.laySsiPlotGeom.addWidget(self.plotWidgetSsiGeom)
            self.laySsiToolGeom.addWidget(NavigationToolbar(self.plotWidgetSsiGeom, self))
            plt.close(_fig)
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setWindowTitle("Scale Error")
            msg.setText(f'Please enter a positive integer')
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()

    # Write results to GUI tables
    def write_to_gui(self, results, table_widget, table_name, name_method):
        headers_labels = []
        if table_name == 'Frequencies':
            table_widget.setColumnCount(self.listIdentifiedPeaks.count())
            for i in range(self.no_channels):
                _temp = f'f_{i + 1}'
                headers_labels.append(_temp)
            table_widget.setHorizontalHeaderLabels(headers_labels)
            table_widget.setRowCount(1)
            data = results['Frequencies']
            if name_method == 'FDD' or name_method == 'SSIcov' or name_method == 'SSIdat':
                for j in range(self.listIdentifiedPeaks.count()):
                    table_widget.setItem(0, j, QTableWidgetItem("{:.6f}".format(data[j])))
            else:
                for j in range(self.listIdentifiedPeaks.count()):
                    table_widget.setItem(0, j, QTableWidgetItem("{:.6f}".format(data[j, 0])))
        elif table_name == 'Damping':
            table_widget.setColumnCount(self.listIdentifiedPeaks.count())
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'd_{i + 1}'
                headers_labels.append(_temp)
            table_widget.setHorizontalHeaderLabels(headers_labels)
            table_widget.setRowCount(1)
            data = results['Damping']
            if name_method == 'FDD' or name_method == 'SSIcov' or name_method == 'SSIdat':
                for j in range(self.listIdentifiedPeaks.count()):
                    table_widget.setItem(0, j, QTableWidgetItem("{:.6f}".format(data[j])))
            else:
                for j in range(self.listIdentifiedPeaks.count()):
                    table_widget.setItem(0, j, QTableWidgetItem("{:.6f}".format(data[j, 0])))
        elif table_name == 'Mode Shapes':
            table_widget.setColumnCount(self.listIdentifiedPeaks.count())
            for i in range(self.listIdentifiedPeaks.count()):
                _temp = f'Mode_{i + 1}'
                headers_labels.append(_temp)
            table_widget.setHorizontalHeaderLabels(headers_labels)
            table_widget.setRowCount(self.no_channels)
            data = results['Mode Shapes'].real
            for i in range(self.no_channels):
                for j in range(self.listIdentifiedPeaks.count()):
                    table_widget.setItem(i, j, QTableWidgetItem("{:.6f}".format(data[i, j])))
                    
    def closeEvent(self, event):
        event.accept()

    # Clicker for clear everything
    def clicker_clear_all(self):
        self.pathFolder = ""
        self.pathLoadFile = ""
        self.inputData = 0
        self.resultsDirectory = ""
        self._figuresFDD = {}
        self._figuresSSI = {}
        self.nodesDict = {}
        self.channelsDict = {}
        self.channelNamesDict = {}
        self.no_channels = 0
        self.nodes = []
        self.connectivity = []
        self.changed_items = []
        self.changed_items_name = []
        self.parameters = Parameters()
        self.fddhelper = 0
        self.my_clicker_helper = 0
        self.modesDict = {}
        self.num_setup = None
        self.buttonLoadData.setEnabled(False)
        self.buttonLoadNodes.setEnabled(False)
        self.buttonLoadConnectivity.setEnabled(False)
        self.buttonCreateGeometry.setEnabled(False)
        self.buttonSubmitChannels.setEnabled(False)
        self.buttonSubmitSetup.setEnabled(False)
        self.buttonAddIdentifiedPeak.setEnabled(False)
        self.buttonDeleteIdentifiedPeak.setEnabled(False)
        self.buttonClearIdentifiedPeaks.setEnabled(False)
        self.buttonSubmitIdentifiedPeaks.setEnabled(False)
        self.buttonFDDSvp.setEnabled(False)
        self.buttonRunFDD.setEnabled(False)
        self.buttonRunSSI.setEnabled(False)
        self.labelCreateFolder.clear()
        self.labelOpenFile.clear()
        self.labelTimeSteps.clear()
        self.labelNumOfChannels.clear()
        self.labelOpenFileNodes.clear()
        self.labelOpenFileConnectivity.clear()
        self.tableNodes.setRowCount(0)
        self.tableDeformedValuesFDD.setRowCount(0)
        self.tableDeformedValuesSSI.setRowCount(0)
        self.listIdentifiedPeaks.clear()
        self.lineEditInsertPeak.clear()
        self.clean_tables(self.tableInputData)
        self.clean_tables(self.tableFreqFdd)
        self.clean_tables(self.tableModeFdd)
        self.clean_tables(self.tableFreqEfdd)
        self.clean_tables(self.tableDampEfdd)
        self.clean_tables(self.tableModeEfdd)
        self.clean_tables(self.tableFreqFsdd)
        self.clean_tables(self.tableDampFsdd)
        self.clean_tables(self.tableModeFsdd)
        self.clean_tables(self.tableFreqSsiCov)
        self.clean_tables(self.tableDampSsiCov)
        self.clean_tables(self.tableModeSsiCov)
        self.clean_tables(self.tableFreqSsiDat)
        self.clean_tables(self.tableDampSsiDat)
        self.clean_tables(self.tableModeSsiDat)
        self.clean_figures(self.layGeometryTool, self.layGeometryPlot, self.plotWidgetGeometry)
        self.clean_figures(self.layFddSvpTool, self.layFddSvpPlot, self.plotWidgetFddSvp)
        self.clean_figures(self.layFddTool, self.layFddPlot, self.plotWidgetFdd)
        self.clean_figures(self.laySsiTool, self.laySsiPlot, self.plotWidgetSsi)
        self.clean_figures(self.layFddToolGeom, self.layFddPlotGeom, self.plotWidgetFddGeom)
        self.clean_figures(self.laySsiToolGeom, self.laySsiPlotGeom, self.plotWidgetSsiGeom)
        self.comboBoxFDDFigures1.clear()
        self.comboBoxFDDFigures2.clear()
        self.comboBoxFDDFiguresGeom1.clear()
        self.comboBoxFDDFiguresGeom2.clear()
        self.comboBoxSSIFigures.clear()
        self.comboBoxSSIFiguresGeom1.clear()
        self.comboBoxSSIFiguresGeom2.clear()
        self.comboBoxChannelNames.clear()
        self.comboBoxChannelAssignments.clear()
        if self.setup_clicks != 1:
            self.tableNodes.itemChanged.disconnect(self.log_change)
        self.setup_clicks = 1
        self.lineEditSamplingFrequency.setText('100')
        self.lineEditDecimationFactor.setText('2')
        self.lineEditTimeShifts.setText('15')
        self.checkBoxDetrend.setChecked(False)
        self.checkBoxDecimation.setChecked(False)
        self.checkBox_OriginalFDD.setChecked(False)
        self.checkBox_FSDD.setChecked(False)
        self.checkBox_EFDD.setChecked(False)
        self.checkBoxDeformedShapeFDD.setChecked(False)
        self.checkBoxValuesOnPlotFDD.setChecked(False)
        self.checkBox_SSI_cov.setChecked(False)
        self.checkBox_SSI_dat.setChecked(False)
        self.checkBoxDeformedShapeSSI.setChecked(False)
        self.checkBoxValuesOnPlotSSI.setChecked(False)
        self.lineEditDimensionScaleFactorXFDD.setText('1')
        self.lineEditDimensionScaleFactorYFDD.setText('1')
        self.lineEditDimensionScaleFactorZFDD.setText('1')
        self.lineEditArrowScaleFactorFDD.setText('1')
        self.lineEditDimensionScaleFactorXSSI.setText('1')
        self.lineEditDimensionScaleFactorYSSI.setText('1')
        self.lineEditDimensionScaleFactorZSSI.setText('1')
        self.lineEditArrowScaleFactorSSI.setText('1')

    @staticmethod
    def remove_widget(lay_tool):
        for i in reversed(range(lay_tool.count())):
            widget_to_remove = lay_tool.itemAt(i).widget()
            # remove it from the layout list
            lay_tool.removeWidget(widget_to_remove)
            # remove it from the gui
            widget_to_remove.setParent(None)

    # Write results to txt
    @staticmethod
    def write_to_txt(results, results_directory, name_method):
        _temp = results_directory + '/' + name_method + '.txt'
        with open(_temp, 'w', encoding='UTF8', newline='') as f:
            f.write(f'Frequencies: ')
            for item in results['Frequencies']:
                f.write("%.6f" % item)
                f.write('  ')
            if name_method != "FDD":
                f.write(f'\nDamping:     ')
                for item in results['Damping']:
                    f.write("%.6f" % item)
                    f.write('  ')
            f.write(f'\n\nMode Shapes: ')
            for item1 in results['Mode Shapes'].real:
                f.write('\n')
                for item in item1:
                    f.write("%.6f" % item)
                    f.write('  ')

    # Clean all tables
    @staticmethod
    def clean_tables(table_widget):
        table_widget.setRowCount(0)
        table_widget.setColumnCount(0)

    # Clean all figures
    @staticmethod
    def clean_figures(lay_tool, lay_plot, plot_widget):
        for i in reversed(range(lay_tool.count())):
            widget_to_remove = lay_tool.itemAt(i).widget()
            lay_tool.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)
        lay_plot.removeWidget(plot_widget)
        _fig = plt.figure(figsize=(8, 6))
        plot_widget = FigureCanvas(_fig)
        lay_plot.addWidget(plot_widget)
        lay_plot.removeWidget(plot_widget)
        plt.clf()

@dataclass
class Parameters:
    sampling_frequency: float = field(default_factory=float)
    decimation_factor: int = field(default_factory=int)
    time_shifts: int = field(default_factory=int)
    identified_peaks: list = field(default_factory=list)
    min_order: int = field(default_factory=int)
    max_order: int = field(default_factory=int)

@dataclass
class Point:
    # xyz[0] --> x, xyz[1] --> y, xyz[2] --> z
    # channels[?] --> Name of the channel in ? direction
    xyz: list[float] = field(default_factory=list)
    channels: list[str] = field(default_factory=list)
    xyz_new: list[float] = field(default_factory=list)


# Run program
def run():
    app = QApplication(sys.argv)
    UIWindow = UI()
    app.exec_()

# TODO: use image metadata to automatically set binning and change region selection
# TODO: add support for an arbitrary number of images, not just K and Rb frames

import sys, os

from PyQt5 import QtWidgets, QtCore
import ctypes

from lib.imageRead import *
from lib.imageProcess import *
from lib.imfitDefaults import *
from lib.gui_helpers import *

import numpy as np

from pymongo import MongoClient, errors
from bson.json_util import loads, dumps

from datetime import datetime


class imfitDue(QtWidgets.QMainWindow):
    def __init__(self, Parent=None):
        super(imfitDue, self).__init__(Parent)

        self.setWindowTitle("ImfitDue: KRb Image Fitting")

        self.regionRb = [0] * 4
        self.regionK = [0] * 4

        self.initializeGui()
        self.createToolbar()
        self.makeConnections()

        self.setupDatabase()

        self.autoloader = Autoloader(self.pf)
        self.autoloader.signalFileArrived.connect(self.loadFile)
        self.autoloader.start()

        self.fitK = None
        self.fitRb = None

        self.mode = DEFAULT_MODE
        self.frame = "OD"

        self.currentFile = None
        self.odK = None
        self.odRb = None

    def setupDatabase(self):
        """
        setupDatabase(self)

        Connects to MongoDB database with configuration specified in ``lib/mongodb.json``.
        """
        try:
            with open("lib/mongodb.json", "r") as f:
                db_config = loads(f.read())
                mongo_url = "mongodb://{}:{}@{}:{}/?authSource=admin".format(
                    db_config["user"],
                    db_config["password"],
                    db_config["address"],
                    db_config["port"],
                )
            self.c = MongoClient(mongo_url, connectTimeoutMS=2000)
            self.c.server_info()
            self.db = self.c["data"]
            self.col = self.db["shots"]
        except Exception as e:
            self.c = None
            self.db = None
            self.col = None
            print("Could not connect to MongoDB: {}".format(e))

    def makeConnections(self):
        self.pf.signalCamChanged.connect(self.camChanged)
        self.pf.loadButton.clicked.connect(self.loadFile)
        self.fo.fitButton.clicked.connect(self.fitCurrent)

        self.figs.plotTools.atomSelectGroup.buttonClicked.connect(self.plotCurrent)
        self.figs.signalFrameChanged.connect(self.frameChanged)

        for i in range(2):
            for j in range(4):
                self.roi.region[i][j].returnPressed.connect(self.currentODCalc)

        self.av.averageButton.clicked.connect(self.averageImages)
        self.fo.uploadButton.clicked.connect(self.process2Origin)
        self.fo.databaseButton.clicked.connect(self.process2Database)

    def camChanged(self, new_cam):
        self.mode = str(new_cam)
        self.autoloader.modeChanged(self.mode)
        self.passCamToROI()

    def frameChanged(self, frame):
        self.frame = str(frame)
        self.plotCurrent()

    def loadFile(self):
        print("loadFile started!")
        path = str(self.pf.filePath.text())
        print("Path read as {}".format(path))
        if os.path.isfile(path):
            print("Path recognized as file!")

            self.autoloader.is_active = False
            self.currentFile = readImage(self.mode, path)

            if "id" in self.currentFile.metadata:
                self.fo.idEdit.setText(self.currentFile.metadata["id"])
                self.fo.idEdit.setStyleSheet(
                    """QLineEdit { background-color: white; }"""
                )
            else:
                self.fo.idEdit.setText("None")
                self.fo.idEdit.setStyleSheet(
                    """QLineEdit { background-color: yellow; }"""
                )

            if self.pf.autoLoad.isChecked():
                t = self.pf.autoLoadFile.text()
                self.pf.autoLoadFile.setText(str(int(t) + 1))

            self.autoloader.is_active = True
            self.currentODCalc()
        else:
            print("File not found!")

    def currentODCalc(self):
        self.autoloader.is_active = False

        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())

        species = IMFIT_MODES[self.mode]["Species"]
        try:
            self.odK = calcOD(self.currentFile, species[0], self.mode, self.regionK)
            self.odRb = calcOD(self.currentFile, species[1], self.mode, self.regionRb)
        except Exception as e:
            print("Could not calculate OD: {}".format(e))

        if self.fo.autoFit.isChecked() and self.frame == "OD":
            self.fitCurrent()
            print("Done fitting current image")
        else:
            self.fitK = None
            self.fitRb = None
            self.plotCurrent()

        self.autoloader.is_active = True
        print("Done calculating current OD")

    def averageImages(self):
        self.autoloader.is_active = False

        self.pf.autoLoad.setChecked(False)
        try:
            x = self.av.getFileNumbers()

            path = IMFIT_MODES[self.mode]["Default Path"]

            firstFile = True
            if x is not None:
                for k in x:
                    print(path.format(k))
                    self.currentFile = readImage(self.mode, path.format(k))
                    if self.currentFile is None:
                        return
                    if firstFile:
                        imageMean = self.currentFile.img
                        firstFile = False
                    else:
                        imageMean += self.currentFile.img

                    self.currentFile.img = imageMean / float(len(x))

                self.currentODCalc()
            self.autoloader.is_active = True
        except Exception as e:
            print("Could not average images: {}".format(e))

    def fitCurrent(self):
        # TODO: Understand what this does and adjust to be more readable
        rbAtom = 1
        kAtom = 0
        TOF = float(self.fo.tof.text())
        pxl = IMFIT_MODES[self.mode]["Pixel Size"]

        fitRbcheckbox = (
            self.fo.fitBothCheckbox.isChecked() and self.fo.fitBothCheckbox.isEnabled()
        ) or "Molecules" not in self.mode

        if self.odK is not None:
            try:
                print(
                    "Fitting K frame: {}".format(
                        str(self.fo.kFitFunction.currentText())
                    )
                )
                self.fitK = fitOD(
                    self.mode,
                    self.odK,
                    str(self.fo.kFitFunction.currentText()),
                    kAtom,
                    TOF,
                    pxl,
                )
                print(processFitResult(self.fitK, self.mode).data_dict)
            except Exception as e:
                self.fitK = None
                print("Could not fit K frame: {}".format(e))
                raise e
        if (self.odRb is not None) and fitRbcheckbox:
            try:
                print(
                    "Fitting Rb frame: {}".format(
                        str(self.fo.rbFitFunction.currentText())
                    )
                )
                self.fitRb = fitOD(
                    self.mode,
                    self.odRb,
                    str(self.fo.rbFitFunction.currentText()),
                    rbAtom,
                    TOF + 6,
                    pxl,
                )
                print(processFitResult(self.fitRb, self.mode).data_dict)
            except Exception as e:
                self.fitRb = None
                print("Could not fit Rb frame: {}".format(e))
                raise e
        else:
            self.fitRb = None
        self.plotCurrent()

        if self.fo.autoUpload.isChecked():
            self.process2Origin()
        if self.fo.autoDatabase.isChecked():
            self.process2Database()
        print("Done fitting and uploading current shot")

    def process2Database(self):
        id = str(self.fo.idEdit.text())
        if id == "None":
            return -1

        result = {"config": IMFIT_MODES[self.mode]}
        if self.fitK is not None:
            KProcess = processFitResult(self.fitK, self.mode)
            species = KProcess.config["Species"][0]
            func = FIT_FUNCTIONS[KProcess.fitObject.fitFunction]
            result[species] = {}
            result[species][func] = KProcess.data_dict
            result[species][func]["region"] = self.regionK

        if self.fitRb is not None:
            RbProcess = processFitResult(self.fitRb, self.mode)
            species = KProcess.config["Species"][1]
            func = FIT_FUNCTIONS[KProcess.fitObject.fitFunction]
            result[species] = {}
            result[species][func] = RbProcess.data_dict
            result[species][func]["region"] = self.regionRb

        try:
            camera_name = self.currentFile.metadata["name"]
        except Exception as e:
            print(
                "Could not set camera name. Assuming name of mode `{}`".format(
                    self.mode
                )
            )
            camera_name = self.mode

        update = [
            {
                "$set": {"images": {camera_name: {"fit": result}}},
            }
        ]
        try:
            self.col.update_one({"_id": id}, update, upsert=True)
        except errors.ConnectionFailure as e:
            print("Could not connect to database: {}\nRetrying connection...".format(e))
            try:
                self.setupDatabase()
                self.col.update_one({"_id": id}, update, upsert=True)
            except Exception as e:
                raise e
        except Exception as e:
            print("Could not upload to database: {}".format(e))
            self.fo.idEdit.setStyleSheet("""QLineEdit { background-color: red; }""")

    def process2Origin(self):
        imagePath = IMFIT_MODES[self.mode]["Image Path"]

        if self.fitK is not None:
            print("Processing K fit result")
            KProcess = processFitResult(self.fitK, self.mode)
            KProcess.data[0] = self.currentFile.fileName + "-" + imagePath
        if self.fitRb is not None:
            print("Processing Rb fit result")
            RbProcess = processFitResult(self.fitRb, self.mode)
            RbProcess.data[0] = self.currentFile.fileName + "-" + imagePath
        if "iXon Molecules" in self.mode:  # Molecule In situ FK
            if (
                self.fo.fitBothCheckbox.isChecked()
                and self.fo.fitBothCheckbox.isEnabled()
            ):
                print("Uploading KRb to Origin")
                upload2Origin(
                    "KRbSpinGauss",
                    self.fitK.fitFunction,
                    [KProcess.data, RbProcess.data],
                )
                print("Done uploading KRb to Origin")
                return 1
            else:
                if "Fermi" in FIT_FUNCTIONS[KProcess.fitObject.fitFunction]:
                    print("Uploading Fermi-Dirac KRb to Origin")
                    upload2Origin("KRb", self.fitK.fitFunction, KProcess.data)
                    print("Done uploading Fermi-Dirac KRb to Origin")
                else:
                    print("Uploading |0,0> KRb to Origin")
                    upload2Origin("KRbFKGauss1", self.fitK.fitFunction, KProcess.data)
                    print("Done uploading |0,0> KRb to Origin")
                return 1
        else:
            if self.fitRb is not None:
                print("Uploading Rb to Origin")
                upload2Origin("Rb", self.fitRb.fitFunction, RbProcess.data)

            print("Uploading K to Origin")
            upload2Origin("K", self.fitK.fitFunction, KProcess.data)
            return 1

    def plotCurrent(self):
        if self.currentFile is None:
            return
        frames = self.currentFile.frames
        species = IMFIT_MODES[self.mode]["Species"]
        if self.figs.plotTools.kSelect.isChecked():

            x = self.odK.xRange0
            y = self.odK.xRange1

            try:
                ch0 = self.fitK.slices.ch0
                ch1 = self.fitK.slices.ch1
                Sx = self.fitK.slices.points0
                Sy = self.fitK.slices.points1
                Fx = self.fitK.slices.fit0
                Fy = self.fitK.slices.fit1

                R = self.fitK.slices.radSlice
                RG = self.fitK.slices.radSliceFitGauss
                RF = self.fitK.slices.radSliceFit
            except Exception as e:
                ch0 = None
                ch1 = None
                Sx = None
                Sy = None
                Fx = None
                Fy = None

                R = None
                RG = None
                RF = None
                self.figs.ax1.cla()
                self.figs.ax2.cla()

                print(e)

            if self.frame == "OD":
                image = self.odK.ODCorrected
            elif self.frame == "Column Density":
                image = self.odK.n
            else:
                image = frames[species[0]][self.frame][y[0] : y[-1], x[0] : x[-1]]
            self.figs.plotUpdate(x, y, image, ch0, ch1)

            if self.frame == "OD":
                if self.fitK is not None:
                    if self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac 2D"):
                        self.figs.plotSliceUpdate(
                            x, [Sx, Fx], np.arange(len(R)), [R, RG, RF]
                        )
                    elif self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac 2D Int"
                    ):
                        self.figs.plotSliceUpdate(x, [Sx, Fx, Fy], y, [Sy])
                    else:
                        self.figs.plotSliceUpdate(x, [Sx, Fx], y, [Sy, Fy])

        if self.figs.plotTools.rbSelect.isChecked():

            x = self.odRb.xRange0
            y = self.odRb.xRange1

            try:
                ch0 = self.fitRb.slices.ch0
                ch1 = self.fitRb.slices.ch1
                Sx = self.fitRb.slices.points0
                Sy = self.fitRb.slices.points1
                Fx = self.fitRb.slices.fit0
                Fy = self.fitRb.slices.fit1
            except Exception as e:
                ch0 = None
                ch1 = None
                Sx = None
                Sy = None
                Fx = None
                Fy = None
                self.figs.ax1.cla()
                self.figs.ax2.cla()
                print(e)

            if self.frame == "OD":
                image = self.odRb.ODCorrected
            elif self.frame == "Column Density":
                image = self.odRb.n
            else:
                image = frames[species[1]][self.frame][y[0] : y[-1], x[0] : x[-1]]
            self.figs.plotUpdate(x, y, image, ch0, ch1)

            if self.frame == "OD":
                self.figs.plotSliceUpdate(x, [Sx, Fx], y, [Sy, Fy])

    def passCamToROI(self):
        self.roi.setDefaultRegion(self.mode)

    def initializeGui(self):

        self.figs = ImageWindows()
        self.fo = fitOptionsWidget(self.figs)
        self.roi = regionWidget()
        self.pf = pathWidget(self.fo, self.figs, self.roi)
        self.av = averageWidget()

        # self.pf.cameraGroup.buttonClicked.connect(self.passCamToROI)

        gb1 = QtWidgets.QGroupBox("File Path")
        gb1.setStyleSheet(self.getStyleSheet("./lib/styles.qss"))
        gb1l = QtWidgets.QVBoxLayout()
        gb1l.addWidget(self.pf)
        gb1.setLayout(gb1l)

        v0 = QtWidgets.QVBoxLayout()
        v0.addStretch(4)
        v0.addWidget(gb1)
        v0.addStretch(1)

        gb2 = QtWidgets.QGroupBox("Region Selection")
        gb2.setStyleSheet(self.getStyleSheet("./lib/styles.qss"))
        gb2l = QtWidgets.QVBoxLayout()
        h0 = QtWidgets.QHBoxLayout()
        h0.addStretch(1)
        h0.addWidget(self.roi)
        h0.addStretch(1)
        gb2l.addLayout(h0)
        gb2.setLayout(gb2l)

        v0.addWidget(gb2)
        v0.addStretch(1)

        gb4 = QtWidgets.QGroupBox("Averaging")
        gb4.setStyleSheet(self.getStyleSheet("./lib/styles.qss"))
        gb4l = QtWidgets.QVBoxLayout()
        gb4l.addWidget(self.av)
        gb4.setLayout(gb4l)

        v0.addWidget(gb4)
        v0.addStretch(1)

        gb3 = QtWidgets.QGroupBox("Fitting Options")
        gb3.setStyleSheet(self.getStyleSheet("./lib/styles.qss"))
        gb3l = QtWidgets.QVBoxLayout()
        gb3l.addWidget(self.fo)
        gb3.setLayout(gb3l)

        v0.addWidget(gb3)

        h = QtWidgets.QHBoxLayout()
        h.addLayout(v0)
        h.addWidget(self.figs)

        self.mainWidget = QtWidgets.QWidget()
        self.mainWidget.setAutoFillBackground(True)
        p = self.mainWidget.palette()
        p.setColor(self.mainWidget.backgroundRole(), QtCore.Qt.white)
        self.mainWidget.setPalette(p)
        self.mainWidget.setLayout(h)

        self.setCentralWidget(self.mainWidget)

    def refreshGui(self):
        # Does not work properly; disabled in menu
        self.initializeGui()

        self.autoloader = Autoloader(self.pf)
        self.autoloader.signalFileArrived.connect(self.loadFile)
        self.autoloader.is_active = True
        self.autoloader.start()

    def saveMainImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(
            self,
            "Save image as",
            "untitled.dat",
            "Data file (*.dat)",
            options=QtWidgets.QFileDialog.DontUseNativeDialog,
        )

        try:
            ok2write = False
            if self.figs.plotTools.kSelect.isChecked():
                if hasattr(self, "odK"):
                    x = self.odK.xRange0
                    y = self.odK.xRange1
                    od = self.odK.ODCorrected
                    ok2write = True
            elif self.figs.plotTools.rbSelect.isChecked():
                if hasattr(self, "odRb"):
                    x = self.odRb.xRange0
                    y = self.odRb.xRange1
                    od = self.odRb.ODCorrected
                    ok2write = True

            if ok2write:
                f = open(xp, "w")
                for i in range(len(y)):
                    for j in range(len(x)):
                        f.write("{0:.3f},".format(od[j, i]))
                    f.write("\n")
                f.close()
            else:
                print("No file to save!")
        except Exception as e:
            print("Could not save image: {}".format(e))

    def saveOSliceImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(
            self,
            "Save image as",
            "untitled.dat",
            "Data file (*.dat)",
            options=QtWidgets.QFileDialog.DontUseNativeDialog,
        )

        try:
            ok2write = False
            if self.figs.plotTools.kSelect.isChecked():
                if hasattr(self, "odK") and hasattr(self.fitK.slices, "points0"):
                    x = self.odK.xRange0
                    y = self.odK.xRange1
                    Sx = self.fitK.slices.points0
                    Fx = self.fitK.slices.fit0
                    Sy = self.fitK.slices.points1
                    Fy = self.fitK.slices.fit1
                    ok2write = True
            elif self.figs.plotTools.rbSelect.isChecked():
                if hasattr(self, "odRb") and hasattr(self.fitRb.slices, "points0"):
                    x = self.odRb.xRange0
                    y = self.odRb.xRange1
                    Sx = self.fitRb.slices.points0
                    Fx = self.fitRb.slices.fit0
                    Sy = self.fitRb.slices.points1
                    Fy = self.fitRb.slices.fit1
                    ok2write = True

            if ok2write:
                f = open(xp, "w")

                xbigger = len(x) >= len(y)
                n = max([len(x), len(y)])
                m = min([len(x), len(y)])

                for k in range(n):
                    if k < m:
                        f.write(
                            "{0:d},{1:.3f},{2:.3f},{3:d},{4:.3f},{5:.3f}\n".format(
                                x[k], Sx[k], Fx[k], y[k], Sy[k], Fy[k]
                            )
                        )
                    elif k > m:
                        if xbigger:
                            f.write(
                                "{0:d},{1:.3f},{2:.3f},,,\n".format(x[k], Sx[k], Fx[k])
                            )
                        elif not xbigger:
                            f.write(
                                ",,,{0:d},{1:.3f},{2:.3f}\n".format(y[k], Sy[k], Fy[k])
                            )

                f.close()
            else:
                print("No file to save!")
        except Exception as e:
            print("Could not save image: {}".format(e))

    def saveRSliceImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(
            self,
            "Save image as",
            "untitled.dat",
            "Data file (*.dat)",
            options=QtWidgets.QFileDialog.DontUseNativeDialog,
        )

        try:
            ok2write = False
            if self.figs.plotTools.kSelect.isChecked():
                if hasattr(self, "odK") and hasattr(self.fitK.slices, "radSlice"):
                    R = self.fitK.slices.radSlice
                    Rf = self.fitK.slices.radSliceFit
                    Rg = self.fitK.slices.radSliceFitGauss
                    ok2write = True
            elif self.figs.plotTools.rbSelect.isChecked():
                if hasattr(self, "odRb") and hasattr(self.fitRb.slices, "radSlice"):
                    R = self.fitRb.slices.radSlice
                    Rf = self.fitRb.slices.radSliceFit
                    Rg = self.fitRb.slices.radSliceFitGauss
                    ok2write = True

            if ok2write:

                f = open(xp, "w")

                for k in range(len(R)):
                    f.write(
                        "{0:.3f},{1:.3f},{2:.3f},{3:.3f}\n".format(
                            k, R[k], Rf[k], Rg[k]
                        )
                    )

                f.close()
            else:
                print("No file to save!")
        except Exception as e:
            print("Could not save image: {}".format(e))

    def loadFromMenu(self):

        d = str(self.pf.filePath.text())
        ext = IMFIT_MODES[self.mode]["Extension Filter"]

        x = QtWidgets.QFileDialog()
        xp = x.getOpenFileName(
            self,
            "Select a file to load",
            filter=ext,
            directory=d,
            options=QtWidgets.QFileDialog.DontUseNativeDialog,
        )
        self.pf.filePath.setText(xp[0])
        self.loadFile()

    def getStyleSheet(self, path):
        f = QtCore.QFile(path)
        f.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text)
        stylesheet = QtCore.QTextStream(f).readAll()
        f.close()

        return stylesheet

    def createToolbar(self):

        exitAction = QtWidgets.QAction("Exit", self)
        exitAction.setShortcut("Ctrl+Q")
        exitAction.triggered.connect(QtWidgets.qApp.quit)

        loadAction = QtWidgets.QAction("Load Image", self)
        loadAction.setShortcut("Ctrl+O")
        loadAction.triggered.connect(self.loadFromMenu)

        # refreshAction = QtWidgets.QAction("Refresh all", self)
        # refreshAction.setShortcut('Ctrl+R')
        # refreshAction.triggered.connect(self.refreshGui)

        saveMain = QtWidgets.QAction("Save Main Image", self)
        saveMain.setShortcut("Ctrl+S")
        saveMain.triggered.connect(self.saveMainImage)
        saveOSlice = QtWidgets.QAction("Save Ortho Slice", self)
        saveOSlice.triggered.connect(self.saveOSliceImage)
        saveRSlice = QtWidgets.QAction("Save Radial Slice", self)
        saveRSlice.triggered.connect(self.saveRSliceImage)

        menubar = self.menuBar()

        fileMenu = menubar.addMenu("File")
        fileMenu.addAction(loadAction)
        # fileMenu.addAction(refreshAction)
        fileMenu.addAction(exitAction)

        saveMenu = menubar.addMenu("Save")
        saveMenu.addAction(saveMain)
        saveMenu.addAction(saveOSlice)
        saveMenu.addAction(saveRSlice)

    def closeEvent(self, event):
        self.autoloader.terminate()


if __name__ == "__main__":
    # The following two lines tell windows that python is only hosting this application
    myappid = "krb.imfitdue"  # arbitrary string
    if os.name == "nt":
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    app = QtWidgets.QApplication(sys.argv)
    w = imfitDue()
    w.setGeometry(100, 100, 1200, 800)

    appico = QtGui.QIcon("main.ico")
    w.setWindowIcon(appico)

    w.show()
    sys.exit(app.exec_())

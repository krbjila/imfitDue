from numpy.ma.core import set_fill_value
import sys

from PyQt5 import QtWidgets, QtCore, QtGui

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt
import numpy as np

from lib.krb_custom_colors import KRbCustomColors
from lib.imfitDefaults import *
from lib.imfitHelpers import *
import os
import datetime


class ImageWindows(QtWidgets.QWidget):
    signalFrameChanged = QtCore.pyqtSignal(str)

    def __init__(self, Parent=None):
        super(ImageWindows, self).__init__(Parent)

        pal = QtGui.QPalette()
        pal.setColor(QtGui.QPalette.Background, QtCore.Qt.white)
        self.setPalette(pal)

        self.frames = {}
        self.species = None
        self.setup()

    def setup(self):
        self.figure = Figure(facecolor="white", tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setContentsMargins(0, 0, 0, 0)
        self.canvas.setFixedSize(1040, 640)

        self.toolbar = NavigationToolbar(self.canvas, self)

        cid = self.canvas.mpl_connect("button_press_event", self.mainGraphClicked)

        ######### Plot Features
        grid = plt.GridSpec(2, 3)
        subplotspec1 = grid.new_subplotspec((0, 0), 2, 2)
        subplotspec2 = grid.new_subplotspec((0, 2), 1, 1)
        subplotspec3 = grid.new_subplotspec((1, 2), 1, 1)

        self.ax0 = self.figure.add_subplot(subplotspec1)
        self.ax1 = self.figure.add_subplot(subplotspec2)
        self.ax2 = self.figure.add_subplot(subplotspec3)

        self.ax0.set_xlabel("X Position")
        self.ax0.set_ylabel("Y Position")
        self.ax1.set_xlabel("X Position")
        self.ax2.set_xlabel("Y Position")

        crossHairColor = [0.5, 0.5, 0.5, 0.5]
        (self.crossHairV,) = self.ax0.plot([], [], color=crossHairColor)
        (self.crossHairH,) = self.ax0.plot([], [], color=crossHairColor)

        ######### Buttons and stuff

        self.plotTools = plotTools()

        self.plotTools.removeCHButton.clicked.connect(self.removeCrossHair)
        self.plotTools.setCHX.returnPressed.connect(self.setCrossHair)
        self.plotTools.setCHY.returnPressed.connect(self.setCrossHair)

        self.plotTools.odSlider.valueChanged.connect(self.plotUpdate)
        self.plotTools.odMinEdit.returnPressed.connect(self.plotUpdate)
        self.plotTools.signalSliderChanged.connect(self.plotUpdate)

        self.plotTools.frameSelect.currentIndexChanged.connect(self.frameChanged)

        ######### Layout Nonsense

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.toolbar)

        vbox2 = QtWidgets.QVBoxLayout()
        vbox2.addWidget(self.plotTools)
        vbox2.addStretch(1)

        box = QtWidgets.QHBoxLayout()
        box.addWidget(self.canvas)
        box.addLayout(vbox2)

        vbox.addLayout(box)

        self.setLayout(vbox)

    def frameChanged(self):
        frame = str(self.plotTools.frameSelect.currentText())
        self.signalFrameChanged.emit(frame)

    def plotUpdate(self, x=None, y=None, image=None, ch0=None, ch1=None, box=None):
        colorMap = KRbCustomColors().whiteJet
        try:
            levelLow = float(self.plotTools.odMinEdit.text())
        except Exception as e:
            levelLow = 0
            print(e)

        if image is not None:
            self.plotTools.setImage(image)

            self.ax0.cla()

            crossHairColor = [0.0, 0.0, 0.0, 1]
            (self.crossHairV,) = self.ax0.plot([], [], color=crossHairColor)
            (self.crossHairH,) = self.ax0.plot([], [], color=crossHairColor)

            try:
                self.mainImage = self.ax0.imshow(
                    image, cmap=colorMap, extent=(min(x), max(x), max(y), min(y))
                )
            except Exception as e:
                print("Could not display image: {}".format(e))
            try:
                self.ax0.set_xlim((min(x), max(x)))
            except Exception as e:
                print("Could not set x limits: {}".format(e))
            try:
                self.ax0.set_ylim((max(y), min(y)))
            except Exception as e:
                print("Could not set y limits: {}".format(e))
            self.setCrossHair()

            if ch0 is not None:
                self.ax0.plot(x, ch0, color=[0.75, 0, 0, 0.75])
            if ch1 is not None:
                self.ax0.plot(ch1, y, color=[0, 0.5, 0, 0.75])
            if box is not None:
                self.ax0.add_patch(
                    Rectangle(
                        (box[0][0], box[1][0]),
                        box[0][1] - box[0][0],
                        box[1][1] - box[1][0],
                        edgecolor="0.5",
                        facecolor="none",
                    )
                )

        try:
            self.mainImage.set_clim(levelLow, self.plotTools.sliderOd())
            self.canvas.draw()
        except Exception as e:
            print("Are you sure you loaded an image?")
            print(e)

    def plotSliceUpdate(self, x, Lx, y, Ly):
        if Lx[0] is not None:
            self.ax1.cla()
            self.ax2.cla()

            for k in range(len(Lx)):
                try:
                    plotStyles = ["ok", "r", "b"]
                    self.ax1.plot(x, Lx[k], plotStyles[k])
                except Exception as e:
                    print("Could not plot x slice: {}".format(e))

            for k in range(len(Ly)):
                try:
                    plotStyles = ["ok", "g", "r"]
                    self.ax2.plot(y, Ly[k], plotStyles[k])
                except Exception as e:
                    print("Could not plot y slice: {}".format(e))

            try:
                self.canvas.draw()
            except Exception as e:
                print("Could not draw plot: {}".format(e))

    def mainGraphClicked(self, event):
        if event.inaxes == self.ax0:

            self.plotTools.setCHX.setText("{0:.1f}".format(event.xdata))
            self.plotTools.setCHY.setText("{0:.1f}".format(event.ydata))

            self.setCrossHair()

    def removeCrossHair(self):
        self.crossHairV.set_xdata([])
        self.crossHairV.set_ydata([])
        self.crossHairH.set_xdata([])
        self.crossHairH.set_ydata([])
        self.canvas.draw()

        self.plotTools.setCHX.setText("")
        self.plotTools.setCHY.setText("")

    def setCrossHair(self):
        xlims = self.ax0.get_xlim()
        ylims = self.ax0.get_ylim()

        try:
            xcoord = float(self.plotTools.setCHX.text())
            ycoord = float(self.plotTools.setCHY.text())
        except:
            return -1

        self.crossHairV.set_xdata([xcoord, xcoord])
        self.crossHairV.set_ydata([ylims[0], ylims[1]])
        self.crossHairH.set_xdata([xlims[0], xlims[1]])
        self.crossHairH.set_ydata([ycoord, ycoord])

        self.ax0.set_xlim(xlims)
        self.ax0.set_ylim(ylims)
        self.canvas.draw()
        return 1


class plotTools(QtWidgets.QWidget):
    signalSliderChanged = QtCore.pyqtSignal(int)

    def __init__(self, Parent=None):
        super(plotTools, self).__init__(Parent)
        self.setFixedHeight(650)
        self.mode = DEFAULT_MODE
        self.image = []
        self.setup()

    def setup(self):
        ODMAXDEFAULT = "2"
        ODMINDEFAULT = "0"
        max_slider = 20
        min_slider = 0

        self.odMaxValidator = QtGui.QDoubleValidator()
        self.odMaxEdit = QtWidgets.QLineEdit(ODMAXDEFAULT)
        self.odMaxEdit.setFixedWidth(60)
        self.odMaxEdit.setValidator(self.odMaxValidator)

        self.odMinValidator = QtGui.QDoubleValidator()
        self.odMinEdit = QtWidgets.QLineEdit(ODMINDEFAULT)
        self.odMinEdit.setFixedWidth(60)
        self.odMinEdit.setValidator(self.odMinValidator)

        self.odMaxValidator.setTop(100)
        self.odMaxValidator.setBottom(float(ODMINDEFAULT))
        self.odMinValidator.setTop(float(ODMAXDEFAULT))
        self.odMinValidator.setBottom(-100)

        self.odMaxEdit.returnPressed.connect(self.updateSlider)
        self.odMinEdit.returnPressed.connect(self.updateSlider)

        self.odSlider = QtWidgets.QSlider(QtCore.Qt.Vertical)
        self.odSlider.setFixedHeight(400)
        self.odSlider.setMaximum(max_slider)
        self.odSlider.setMinimum(min_slider)
        self.odSlider.setValue(10)
        self.odSlider.setTickPosition(2)
        self.odSlider.setTickInterval(1)
        self.odSlider.setPageStep(2)

        self.removeCHButton = QtWidgets.QPushButton("Remove\nCrosshair")

        self.setCHX = QtWidgets.QLineEdit("0")
        self.setCHY = QtWidgets.QLineEdit("0")

        self.atomSelectGroup = QtWidgets.QButtonGroup()
        self.rbSelect = QtWidgets.QRadioButton("Rb")
        self.kSelect = QtWidgets.QRadioButton("K")
        self.kSelect.setChecked(True)
        self.atomSelectGroup.addButton(self.rbSelect)
        self.atomSelectGroup.addButton(self.kSelect)

        self.frameSelect = QtWidgets.QComboBox()
        self.frameSelect.addItem("OD")
        self.frameSelect.addItem("Shadow")
        self.frameSelect.addItem("Light")
        self.frameSelect.addItem("Dark")
        self.frameSelect.addItem("Column Density")
        self.frameSelect.setCurrentIndex(0)

        self.autoscaler = QtWidgets.QPushButton("Autoscale me!")
        self.autoscaler.clicked.connect(self.autoscaleSlider)

        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(self.setCHX)
        hbox.addWidget(self.setCHY)

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.odMaxEdit, QtCore.Qt.AlignCenter)
        vbox.addWidget(self.odSlider)
        vbox.addWidget(self.odMinEdit)
        vbox.addWidget(self.autoscaler)
        vbox.addWidget(self.removeCHButton)
        vbox.addLayout(hbox)
        vbox.addWidget(QtWidgets.QLabel("Crosshair     (X,Y)"))

        vbox.addStretch(1)
        vbox.addWidget(self.kSelect)
        vbox.addWidget(self.rbSelect)

        vbox.addWidget(self.frameSelect)

        self.setLayout(vbox)

    def sliderOd(self):
        return float(self.odMinEdit.text()) + float(self.odMaxEdit.text()) * (
            self.odSlider.value() - self.odSlider.minimum()
        ) / (self.odSlider.maximum() - self.odSlider.minimum())

    def updateSlider(self):
        odmax = self.odMaxEdit.text()
        odmin = self.odMinEdit.text()
        sv = self.sliderOd()

        try:
            odmax = float(odmax)
            odmin = float(odmin)
        except:
            return -1

        self.odMaxValidator.setBottom(odmin)
        self.odMinValidator.setTop(odmax)

        if odmax < sv:
            self.odSlider.setValue(self.odSlider.maximum())

        if odmin > sv:
            self.odSlider.setValue(self.odSlider.minimum())

        self.signalSliderChanged.emit(self.odSlider.value())

    def setImage(self, image):
        self.image = image

    def autoscaleSlider(self):
        if len(self.image):
            self.odMaxEdit.setValidator(None)
            self.odMinEdit.setValidator(None)
            (low, high) = np.percentile(self.image, [AUTOSCALE_MIN, 100])
            self.odMaxEdit.setText("{:0.2f}".format(high * AUTOSCALE_HEADROOM))
            self.odMinEdit.setText("{:0.2f}".format(low))
            self.odSlider.setValue(self.odSlider.maximum())
            self.odMaxEdit.setValidator(self.odMaxValidator)
            self.odMinEdit.setValidator(self.odMinValidator)
            self.updateSlider()


class regionWidget(QtWidgets.QWidget):
    def __init__(self, Parent=None):
        super(regionWidget, self).__init__(Parent)
        self.setup()

    def setup(self):

        topLabels = ["XC", "YC", "CrX", "CrY"]
        # sideLabels = ATOM_NAMES
        sideLabels = IMFIT_MODES[DEFAULT_MODE]["Species"]

        font = QtGui.QFont()
        font.setBold(True)
        font.setPointSize(12)

        self.grid = QtWidgets.QGridLayout()
        self.atom_labels = []
        for k in range(4):
            x = QtWidgets.QLabel(topLabels[k])
            x.setFont(font)
            self.grid.addWidget(x, 0, k + 1, 1, 1)
        for k in range(2):
            x = QtWidgets.QLabel(sideLabels[k])
            self.atom_labels.append(x)
            x.setFont(font)
            self.grid.addWidget(x, k + 1, 0, 1, 1)

        self.region = [[0] * 4, [0] * 4]

        for i in range(2):
            for j in range(4):
                self.region[i][j] = QtWidgets.QLineEdit(
                    str(IMFIT_MODES[DEFAULT_MODE]["Default Region"][i][j])
                )
                self.region[i][j].setValidator(QtGui.QIntValidator())
                self.region[i][j].setFixedWidth(50)
                self.grid.addWidget(self.region[i][j], i + 1, j + 1, 1, 1)

        self.setLayout(self.grid)

    def setDefaultRegion(self, mode):
        try:
            region = IMFIT_MODES[mode]["Default Region"]
        except AttributeError:
            print("regionWidget.setDefaultRegion error: requested mode doesn't exist")
            return -1

        for i in range(2):
            for j in range(4):
                self.region[i][j].setText(str(region[i][j]))
        return 0


class pathWidget(QtWidgets.QWidget):
    signalCamChanged = QtCore.pyqtSignal(str)

    def __init__(self, fitOptions, imageWindows, roi, Parent=None):
        super(pathWidget, self).__init__(Parent)
        self.fitOptions = fitOptions
        self.imageWindows = imageWindows
        self.roi = roi
        self.mode = DEFAULT_MODE
        self.setup()

    def setup(self):

        self.filePath = QtWidgets.QLineEdit(IMFIT_MODES[self.mode]["Default Path"])
        self.filePath.setFixedWidth(300)
        self.browseButton = QtWidgets.QPushButton("Browse")
        self.browseButton.clicked.connect(self.browseFile)
        self.loadButton = QtWidgets.QPushButton("Load")
        self.autoLoadFile = QtWidgets.QLineEdit("0")
        self.autoLoadFile.setValidator(QtGui.QIntValidator())
        self.autoLoadFile.setFixedWidth(60)
        self.autoLoad = QtWidgets.QCheckBox("AutoLoad")

        self.cameraGroup = QtWidgets.QComboBox()
        self.cameraGroup.addItems(IMFIT_MODES.keys())
        self.cameraGroup.setCurrentIndex(0)
        self.cameraGroup.currentIndexChanged.connect(self.camChanged)

        h0 = QtWidgets.QHBoxLayout()
        h1 = QtWidgets.QHBoxLayout()
        h2 = QtWidgets.QHBoxLayout()
        v = QtWidgets.QVBoxLayout()

        # h0.addStretch(1)
        # h0.addWidget(self.pi)
        # h0.addWidget(self.ixon)
        # h0.addWidget(self.ixon_gsm)
        # h0.addWidget(self.ixonv)
        h0.addStretch(1)
        h0.addWidget(QtWidgets.QLabel("Imaging/Analysis mode: "))
        h0.addWidget(self.cameraGroup)

        h1.addWidget(self.filePath)
        h1.addWidget(self.browseButton)
        h1.addWidget(self.loadButton)

        h2.addStretch(1)
        h2.addWidget(QtWidgets.QLabel("File #:"))
        h2.addWidget(self.autoLoadFile)
        h2.addWidget(self.autoLoad)

        v.addLayout(h0)
        v.addLayout(h1)
        v.addLayout(h2)

        self.setLayout(v)
        self.camChanged()

    def browseFile(self):

        d = self.filePath.text()
        if not os.path.isdir(d):
            d = IMFIT_MODES[self.mode]["Default Path"]

        ext = IMFIT_MODES[self.mode]["Extension Filter"]

        x = QtWidgets.QFileDialog()
        (path, _) = x.getOpenFileName(
            self,
            "Select a file to load",
            filter=ext,
            directory=d,
            options=QtWidgets.QFileDialog.DontUseNativeDialog,
        )
        self.filePath.setText(path)

    def camChanged(self):
        self.autoLoad.setChecked(False)

        self.mode = str(self.cameraGroup.currentText())
        self.signalCamChanged.emit(self.mode)  # sets self.mode in top-level class

        self.fitOptions.rbFitFunction.clear()
        self.fitOptions.rbFitFunction.addItems(IMFIT_MODES[self.mode]["Fit Functions"])
        self.fitOptions.kFitFunction.clear()
        self.fitOptions.kFitFunction.addItems(IMFIT_MODES[self.mode]["Fit Functions"])
        if IMFIT_MODES[self.mode]["Enforce same fit for both"]:
            kfit = self.fitOptions.kFitFunction.currentIndex()
            self.fitOptions.rbFitFunction.setCurrentIndex(kfit)
            self.fitOptions.rbFitFunction.setEnabled(False)
        else:
            self.fitOptions.rbFitFunction.setEnabled(True)
        if "Allow fit both states" in IMFIT_MODES[self.mode]:
            self.fitOptions.fitBothCheckbox.setEnabled(True)
            if not self.fitOptions.fitBothCheckbox.isChecked():
                self.imageWindows.plotTools.rbSelect.setEnabled(False)
                self.imageWindows.plotTools.kSelect.setChecked(True)
            else:
                self.imageWindows.plotTools.rbSelect.setEnabled(True)
        else:
            self.fitOptions.fitBothCheckbox.setEnabled(False)
            self.imageWindows.plotTools.rbSelect.setEnabled(True)
        self.fitOptions.KLabel.setText(
            "Fit {} to".format(IMFIT_MODES[self.mode]["Species"][0])
        )
        self.fitOptions.RbLabel.setText(
            "Fit {} to".format(IMFIT_MODES[self.mode]["Species"][1])
        )
        self.imageWindows.plotTools.kSelect.setText(
            IMFIT_MODES[self.mode]["Species"][0]
        )
        self.imageWindows.plotTools.rbSelect.setText(
            IMFIT_MODES[self.mode]["Species"][1]
        )
        self.roi.atom_labels[0].setText(IMFIT_MODES[self.mode]["Species"][0])
        self.roi.atom_labels[1].setText(IMFIT_MODES[self.mode]["Species"][1])

        d = IMFIT_MODES[self.mode]["Default Path"]
        self.filePath.setText(d)
        if os.path.isdir(d):
            n = getLastFile(d, IMFIT_MODES[self.mode]["Extension Filter"])
            self.autoLoadFile.setText(str(n))


class fitOptionsWidget(QtWidgets.QWidget):
    def __init__(self, plots, Parent=None):
        super(fitOptionsWidget, self).__init__(Parent)
        self.plots = plots
        self.setup()

    def updateRbFit(self):
        if not self.rbFitFunction.isEnabled():
            self.rbFitFunction.setCurrentIndex(self.kFitFunction.currentIndex())

    def updateRbEnabled(self):
        if self.fitBothCheckbox.isEnabled() and not self.fitBothCheckbox.isChecked():
            self.plots.plotTools.rbSelect.setEnabled(False)
            self.plots.plotTools.kSelect.setChecked(True)
        else:
            self.plots.plotTools.rbSelect.setEnabled(True)

    def setup(self):
        self.rbFitFunction = QtWidgets.QComboBox()
        self.rbFitFunction.addItems(IMFIT_MODES[DEFAULT_MODE]["Fit Functions"])
        self.kFitFunction = QtWidgets.QComboBox()
        self.kFitFunction.addItems(IMFIT_MODES[DEFAULT_MODE]["Fit Functions"])
        self.kFitFunction.currentIndexChanged.connect(self.updateRbFit)

        # self.imagePath = QtWidgets.QComboBox()
        # self.imagePath.addItems(IMAGING_PATHS)

        self.autoFit = QtWidgets.QCheckBox("AutoFit")
        self.autoUpload = QtWidgets.QCheckBox("Auto Origin")
        self.fitBothCheckbox = QtWidgets.QCheckBox("Fit both states?")
        self.fitBothCheckbox.setEnabled(
            "Allow fit both states" in IMFIT_MODES[DEFAULT_MODE]
        )
        self.fitBothCheckbox.stateChanged.connect(self.updateRbEnabled)

        self.fitButton = QtWidgets.QPushButton("Fit")
        self.uploadButton = QtWidgets.QPushButton("Upload to Origin")

        self.autoDatabase = QtWidgets.QCheckBox("Auto Database")
        self.idEdit = QtWidgets.QLineEdit("None")
        self.idEdit.setDisabled(True)
        self.databaseButton = QtWidgets.QPushButton("Upload to Database")

        self.tof = QtWidgets.QLineEdit("6")
        tof_validator = QtGui.QDoubleValidator()
        tof_validator.setBottom(0)
        self.tof.setValidator(tof_validator)
        self.tof.setFixedWidth(30)

        #### Layout Stuff

        h0 = QtWidgets.QHBoxLayout()
        h0.addStretch(1)
        self.RbLabel = QtWidgets.QLabel("Fit Rb to:")
        h0.addWidget(self.RbLabel)
        h0.addWidget(self.rbFitFunction)

        h1 = QtWidgets.QHBoxLayout()
        h1.addWidget(QtWidgets.QLabel("K TOF (ms):"))
        h1.addWidget(self.tof)
        h1.addStretch(1)
        self.KLabel = QtWidgets.QLabel("Fit K to:")
        h1.addWidget(self.KLabel)
        h1.addWidget(self.kFitFunction)

        h1b = QtWidgets.QHBoxLayout()
        h1b.addStretch(1)
        self.idLabel = QtWidgets.QLabel("ID:")
        h1b.addWidget(self.idLabel)
        h1b.addWidget(self.idEdit)

        h2 = QtWidgets.QHBoxLayout()
        h2.addWidget(self.fitButton)
        h2.addWidget(self.uploadButton)
        h2.addWidget(self.databaseButton)

        v0 = QtWidgets.QVBoxLayout()
        v0.addLayout(h0)
        v0.addLayout(h1)
        v0.addLayout(h1b)
        v0.addLayout(h2)

        v1 = QtWidgets.QGridLayout()
        v1.addWidget(self.autoFit)
        v1.addWidget(self.autoUpload)
        v1.addWidget(self.autoDatabase)
        v1.addWidget(self.fitBothCheckbox)

        h = QtWidgets.QHBoxLayout()
        h.addLayout(v0)
        h.addLayout(v1)

        self.setLayout(h)


class averageWidget(QtWidgets.QWidget):
    def __init__(self, Parent=None):
        super(averageWidget, self).__init__(Parent)
        self.setup()

    def setup(self):

        self.averageEdit = QtWidgets.QLineEdit()
        self.averageButton = QtWidgets.QPushButton("Average")

        box = QtWidgets.QHBoxLayout()
        box.addWidget(QtWidgets.QLabel("Images to Average: "))
        box.addWidget(self.averageEdit)
        box.addWidget(self.averageButton)

        self.setLayout(box)

    def getFileNumbers(self):
        x = self.averageEdit.text()
        if x == "":
            return None
        else:
            return getImagesFromRange(x)


class Autoloader(QtCore.QThread):
    # TODO: Redo this to support autoloading differently based on self.modes (for creating class, not this one)

    signalFileArrived = QtCore.pyqtSignal()

    def __init__(self, mainPF, Parent=None):
        super(Autoloader, self).__init__(Parent)

        self.mainPF = mainPF
        self.startDate = datetime.datetime.now().strftime("%d")
        self.autoloadState = True

        self.mode = DEFAULT_MODE
        self.config = IMFIT_MODES

        self.is_active = True

    def wait_a_while(self):
        self.msleep(1000)

    def modeChanged(self, mode):
        self.mode = mode

    def run(self):
        from os import path

        inner_wait = 1000
        outer_wait = 1000

        while True:
            if self.mainPF.autoLoad.isChecked() and self.is_active:
                if self.startDate != datetime.datetime.now().strftime("%d"):
                    import lib.imfitDefaults as imfitDefaults  # Force reload date in path

                    self.config = imfitDefaults.IMFIT_MODES
                    self.startDate = datetime.datetime.now().strftime("%d")

                try:
                    nextFile = int(self.mainPF.autoLoadFile.text())
                except ValueError:
                    print("Something other than an integer is in the autoload box!")

                # print("Checking file {}".format(nextFile))

                nextPath = self.config[self.mode]["Default Path"] + self.config[
                    self.mode
                ]["Default Suffix"].format(nextFile)
                if path.isfile(nextPath):
                    self.msleep(inner_wait)
                    self.mainPF.filePath.setText(nextPath)
                    self.signalFileArrived.emit()
            self.msleep(outer_wait)


def main():
    app = QtWidgets.QApplication(sys.argv)

    w = averageWidget()

    w.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

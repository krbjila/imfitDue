import sys

from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar 
from matplotlib.figure import Figure 

import matplotlib.pyplot as plt
import numpy as np

from krb_custom_colors import KRbCustomColors
from imfitDefaults import *
from imfitHelpers import *
import os


class ImageWindows(QtGui.QWidget):
    def __init__(self, Parent=None):
        super(ImageWindows,self).__init__(Parent)
        # self.setGeometry(10,10,600,400)
        # self.setGeometry(10,10,1000,600)
        
        pal = QtGui.QPalette()
        pal.setColor(QtGui.QPalette.Background, QtCore.Qt.white)
        self.setPalette(pal)
        
        
        self.setup()


    def setup(self):
        self.figure = Figure(facecolor='white',tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setContentsMargins(0,0,0,0)
        self.canvas.setFixedSize(1040,640)

        self.toolbar = NavigationToolbar(self.canvas, self)

        cid = self.canvas.mpl_connect('button_press_event', self.mainGraphClicked)

        ######### Plot Features 
        grid = plt.GridSpec(2,3)
        subplotspec1 = grid.new_subplotspec((0,0),2,2)
        subplotspec2 = grid.new_subplotspec((0,2),1,1)
        subplotspec3 = grid.new_subplotspec((1,2),1,1)

        self.ax0 = self.figure.add_subplot(subplotspec1)
        self.ax1 = self.figure.add_subplot(subplotspec2)
        self.ax2 = self.figure.add_subplot(subplotspec3)

        self.ax0.set_xlabel('X Position')
        self.ax0.set_ylabel('Y Position')
        self.ax1.set_xlabel('X Position')
        self.ax2.set_xlabel('Y Position')


        crossHairColor = [0.5, 0.5, 0.5,0.5]
        self.crossHairV, = self.ax0.plot([],[],color=crossHairColor)
        self.crossHairH, = self.ax0.plot([],[],color=crossHairColor)

        ######### Buttons and stuff

        self.plotTools = plotTools()

        self.plotTools.removeCHButton.clicked.connect(self.removeCrossHair)
        self.plotTools.setCHX.returnPressed.connect(self.setCrossHair)
        self.plotTools.setCHY.returnPressed.connect(self.setCrossHair)

        self.plotTools.odSlider.valueChanged.connect(self.plotUpdate)
        self.plotTools.odMinEdit.returnPressed.connect(self.plotUpdate)


        ######### Layout Nonsense

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.toolbar)

        vbox2 = QtGui.QVBoxLayout()
        vbox2.addWidget(self.plotTools)
        vbox2.addStretch(1)


        box = QtGui.QHBoxLayout()
        box.addWidget(self.canvas)
        box.addLayout(vbox2)
        
        vbox.addLayout(box)

        self.setLayout(vbox)



    def plotUpdate(self, x=None, y=None, image=None, ch0=None, ch1=None):



        colorMap = KRbCustomColors().whiteJet
        levelLow = float(self.plotTools.odMinEdit.text())
        levelHigh = float(self.plotTools.odMaxEdit.text())

        if image is not None:
            self.ax0.cla()

            crossHairColor = [0.5, 0.5, 0.5,0.5]
            self.crossHairV, = self.ax0.plot([],[],color=crossHairColor)
            self.crossHairH, = self.ax0.plot([],[],color=crossHairColor)
            
            self.mainImage = self.ax0.pcolormesh(x,y,image,cmap=colorMap)
            self.setCrossHair()            
            
            if ch0 is not None:
                self.ax0.plot(x,ch0,color=[1,0,0,0.25])
            if ch1 is not None:
                self.ax0.plot(ch1,y,color=[0,1,0,0.25])
        
        try:
            self.mainImage.set_clim(levelLow, self.plotTools.odSlider.value()/10.0)
        except:
            print('Are you sure you loaded an image?')

        self.canvas.draw()

    def plotSliceUpdate(self, x, Lx, y, Ly):
        if Lx[0] is not None:
            self.ax1.cla()
            self.ax2.cla()

            plotStyles = ['ok', 'r', 'b']

            for k in range(len(Lx)):
                self.ax1.plot(x,Lx[k],plotStyles[k])
            
            for k in range(len(Ly)):
                 self.ax2.plot(y,Ly[k],plotStyles[k])

            self.canvas.draw()


    def mainGraphClicked(self, event):
        if event.inaxes == self.ax0:

            self.plotTools.setCHX.setText('{0:.1f}'.format(event.xdata))
            self.plotTools.setCHY.setText('{0:.1f}'.format(event.ydata))

            self.setCrossHair()
#            self.plotUpdate()

    def removeCrossHair(self):
        self.crossHairV.set_xdata([])
        self.crossHairV.set_ydata([])
        self.crossHairH.set_xdata([])
        self.crossHairH.set_ydata([])
        self.canvas.draw()

        self.plotTools.setCHX.setText('')
        self.plotTools.setCHY.setText('')

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


class plotTools(QtGui.QWidget):
    def __init__(self,Parent=None):
        super(plotTools,self).__init__(Parent)
        self.setup()
        self.setFixedHeight(650)

    def setup(self):

        ODMAXDEFAULT = '2'
        ODMINDEFAULT = '0'


        self.odMaxEdit = QtGui.QLineEdit(ODMAXDEFAULT)
        self.odMaxEdit.setFixedWidth(60)

        self.odMinEdit = QtGui.QLineEdit(ODMINDEFAULT)
        self.odMinEdit.setFixedWidth(60)

        self.odMaxEdit.returnPressed.connect(self.updateSlider)
        self.odMinEdit.returnPressed.connect(self.updateSlider)
        
        self.odSlider = QtGui.QSlider(QtCore.Qt.Vertical)
        self.odSlider.setFixedHeight(400)
        self.odSlider.setMaximum(float(ODMAXDEFAULT)*10.0)
        self.odSlider.setMinimum(float(ODMINDEFAULT)*10.0+0.5)
        self.odSlider.setValue(10)
        self.odSlider.setTickPosition(2)
        self.odSlider.setTickInterval(5)
        self.odSlider.setPageStep(2)

        self.removeCHButton = QtGui.QPushButton('Remove\nCrosshair')

        self.setCHX = QtGui.QLineEdit('0')
        self.setCHY = QtGui.QLineEdit('0')

        self.atomSelectGroup = QtGui.QButtonGroup()
        self.rbSelect = QtGui.QRadioButton('Rb')
        self.kSelect = QtGui.QRadioButton('K')
        self.kSelect.setChecked(True)
        self.atomSelectGroup.addButton(self.rbSelect)
        self.atomSelectGroup.addButton(self.kSelect)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.setCHX)
        hbox.addWidget(self.setCHY)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.odMaxEdit,QtCore.Qt.AlignCenter)
        vbox.addWidget(self.odSlider)
        vbox.addWidget(self.odMinEdit)
        vbox.addWidget(self.removeCHButton)
        vbox.addLayout(hbox)
        vbox.addWidget(QtGui.QLabel('Crosshair     (X,Y)'))

        vbox.addStretch(1)
        vbox.addWidget(self.kSelect)
        vbox.addWidget(self.rbSelect)


        #self.addLayout(vbox)
        self.setLayout(vbox)

    def updateSlider(self):
        
        odmax = self.odMaxEdit.text()
        odmin = self.odMinEdit.text()
        sv = self.odSlider.value()

        try:
            odmax = float(odmax)*10
            odmin = float(odmin)*10
        except:

            return -1

        if odmax < sv:
            self.odSlider.setValue(odmax)
        self.odSlider.setMaximum(odmax)

        if odmin > sv:
            self.odSlider.setValue(odmin)
        self.odSlider.setMinimum(odmin+0.5)
            

class regionWidget(QtGui.QWidget):
    def __init__(self, Parent=None):
        super(regionWidget, self).__init__(Parent)
        self.setup()

    def setup(self):

        topLabels = ['XC', 'YC', 'CrX', 'CrY']
        sideLabels = ATOM_NAMES



        font = QtGui.QFont()
        font.setBold(True)
        font.setPointSize(12)

        grid = QtGui.QGridLayout()
        for k in range(4):
            x = QtGui.QLabel(topLabels[k])
            x.setFont(font)
            grid.addWidget(x,0,k+1,1,1)
        for k in range(2):
            x = QtGui.QLabel(ATOM_NAMES[k])
            x.setFont(font)
            grid.addWidget(x,k+1,0,1,1)
        
        self.region = [[0]*4,[0]*4]
        for i in range(2):
            for j in range(4):
                self.region[i][j] = QtGui.QLineEdit(str(DEFAULT_REGION[i][j]))
                self.region[i][j].setFixedWidth(50)
                grid.addWidget(self.region[i][j],i+1,j+1,1,1)

        self.setLayout(grid)

class pathWidget(QtGui.QWidget):
    def __init__(self,Parent=None):
        super(pathWidget,self).__init__(Parent)
        self.setup()

    def setup(self):


        self.filePath = QtGui.QLineEdit(DEFAULT_PATH)
        self.filePath.setFixedWidth(300)
        self.browseButton = QtGui.QPushButton('Browse')
        self.browseButton.clicked.connect(self.browseFile)
        self.loadButton = QtGui.QPushButton('Load')
        self.autoLoadFile = QtGui.QLineEdit('0')
        self.autoLoadFile.setFixedWidth(60)
        self.autoLoad = QtGui.QCheckBox('AutoLoad')


        self.cameraGroup = QtGui.QButtonGroup()
        self.pi = QtGui.QRadioButton('Princeton Instruments')
        self.pi.setChecked(True)
        self.ixon = QtGui.QRadioButton('Andor iXon 888')
        self.cameraGroup.addButton(self.pi, 0)
        self.cameraGroup.addButton(self.ixon, 1)
        self.cameraGroup.buttonClicked.connect(self.camChanged)        
        
        h0 = QtGui.QHBoxLayout()
        h1  = QtGui.QHBoxLayout()
        h2  = QtGui.QHBoxLayout()
        v  = QtGui.QVBoxLayout()

        h0.addStretch(1)
        h0.addWidget(QtGui.QLabel('Camera: '))
        h0.addWidget(self.pi)
        h0.addWidget(self.ixon)


        h1.addWidget(self.filePath)
        h1.addWidget(self.browseButton)
        h1.addWidget(self.loadButton)

        h2.addStretch(1)
        h2.addWidget(QtGui.QLabel('File #:'))
        h2.addWidget(self.autoLoadFile)
        h2.addWidget(self.autoLoad)

        v.addLayout(h0)
        v.addLayout(h1)
        v.addLayout(h2)
        
        
        self.setLayout(v)

    def browseFile(self):
        
        d = self.filePath.text()
        if not os.path.isdir(d):
            d = DEFAULT_PATH        
        
        
        if self.cameraGroup.checkedId():
            ext = '*.csv'
        else:
            ext = '*.spe'


        x = QtGui.QFileDialog()
        path = x.getOpenFileName(self,'Select a file to load' ,filter=ext, directory=d, options=QtGui.QFileDialog.DontUseNativeDialog)
        self.filePath.setText(path)


    def camChanged(self):
        self.autoLoad.setChecked(False)
        if self.cameraGroup.checkedId():
            d = DEFAULT_PATH_IXON
            self.filePath.setText(d)
            if os.path.isdir(d):
                n = getLastFile(d)
                self.autoLoadFile.setText(str(n))
        else:
            d = DEFAULT_PATH_PI
            self.filePath.setText(d)
            if os.path.isdir(d):
                n = getLastFile(d)
                self.autoLoadFile.setText(str(n))

class fitOptionsWidget(QtGui.QWidget):
    def __init__(self,Parent=None):
        super(fitOptionsWidget,self).__init__(Parent)
        self.setup()

    def setup(self):


        self.rbFitFunction = QtGui.QComboBox()
        self.rbFitFunction.addItems(FIT_FUNCTIONS)
        self.kFitFunction = QtGui.QComboBox()
        self.kFitFunction.addItems(FIT_FUNCTIONS)

        self.imagePath = QtGui.QComboBox()
        self.imagePath.addItems(IMAGING_PATHS)

        self.autoFit = QtGui.QCheckBox('AutoFit')
        self.autoUpload = QtGui.QCheckBox('Auto Origin')
        self.moleculeBook = QtGui.QCheckBox('Molecules?')

        self.fitButton = QtGui.QPushButton('Fit')
        self.uploadButton = QtGui.QPushButton('Upload to Origin')

        
        #### Layout Stuff

        h0 = QtGui.QHBoxLayout()
        h0.addWidget(QtGui.QLabel("Imageing Path: "))
        h0.addWidget(self.imagePath)
        h0.addStretch(1)
        h0.addWidget(QtGui.QLabel('Fit Rb to:'))
        h0.addWidget(self.rbFitFunction)

        h1 = QtGui.QHBoxLayout()
        h1.addStretch(1)
        h1.addWidget(QtGui.QLabel('Fit K to:'))
        h1.addWidget(self.kFitFunction)

        h2 = QtGui.QHBoxLayout()
        h2.addWidget(self.fitButton)
        h2.addWidget(self.uploadButton)

        v0 = QtGui.QVBoxLayout()
        v0.addLayout(h0)
        v0.addLayout(h1)
        v0.addLayout(h2)

        v1 = QtGui.QGridLayout()
        v1.addWidget(self.autoFit)
        v1.addWidget(self.autoUpload)
        v1.addWidget(self.moleculeBook)

        h = QtGui.QHBoxLayout()
        h.addLayout(v0)
        h.addLayout(v1)

        self.setLayout(h)


class averageWidget(QtGui.QWidget):
    def __init__(self, Parent=None):
        super(averageWidget,self).__init__(Parent) 
        self.setup()

    def setup(self):

        self.averageEdit = QtGui.QLineEdit() 
        self.averageButton = QtGui.QPushButton("Average")

        box = QtGui.QHBoxLayout()
        box.addWidget(QtGui.QLabel("Images to Average: "))
        box.addWidget(self.averageEdit)
        box.addWidget(self.averageButton)

        self.setLayout(box)



    def getFileNumbers(self): 
        x = self.averageEdit.text()
        if x == '':
            return None
        else:
            return getImagesFromRange(x)



class autoloader(QtCore.QThread):


    def __init__(self, mainPF, Parent=None):
        super(autoloader, self).__init__(Parent)

        self.mainPF = mainPF

        self.autoloadState = True
        self.camera = None

    def run(self):
        from os import path 

        while True:
            if self.mainPF.autoLoad.isChecked():

                camera = self.mainPF.cameraGroup.checkedId()

                nextFile = self.mainPF.autoLoadFile.text()

                fileGood = 0
                try:
                    int(nextFile)
                    fileGood = 1
                except:
                    pass

                if camera == 0 and fileGood == 1:
                    nextPath = DEFAULT_PATH_PI + "pi_" + nextFile + ".spe"
                    if path.isfile(nextPath):
                    	self.mainPF.filePath.setText(nextPath)
                        self.emit(QtCore.SIGNAL('fileArrived'))

                elif camera == 1 and fileGood == 1:
                    nextPath = DEFAULT_PATH_IXON + "ixon_" + nextFile + ".csv"
                    if path.isfile(nextPath):
                        self.mainPF.filePath.setText(nextPath)
                        self.emit(QtCore.SIGNAL('fileArrived'))

            self.sleep(1)

def main():
    app = QtGui.QApplication(sys.argv)

    w = averageWidget()

    w.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

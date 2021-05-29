from lib.imfitDefaults import IMFIT_MODES
import sys, os
sys.path.append('./lib')

from PyQt5 import QtWidgets, QtCore
import ctypes

from imageRead import *
from imageProcess import *
from imfitDefaults import *
from gui_helpers import *

import numpy as np

class imfitDue(QtWidgets.QMainWindow):
    def __init__(self,Parent=None):
        super(imfitDue, self).__init__(Parent)

        self.setWindowTitle('ImfitDue: KRb Image Fitting')

        self.regionRb = [0]*4
        self.regionK = [0]*4


        self.initializeGui()
        self.createToolbar()
        self.makeConnections()

        self.autoloader = Autoloader(self.pf)
        self.autoloader.signalFileArrived.connect(self.loadFile)
        # self.autoloader.connect(self.autoloader, QtCore.SIGNAL('fileArrived'), self.loadFile)
        self.autoloader.start()

        self.fitK = None
        self.fitRb = None

        self.mode = DEFAULT_MODE
        self.frame = 'OD'

        self.currentFile = None
        self.odK = None
        self.odRb = None

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

            # TODO: Implement readImage
            self.currentFile = readImage(self.mode, path)
            
            if self.pf.autoLoad.isChecked():
                t = self.pf.autoLoadFile.text()
                self.pf.autoLoadFile.setText(str(int(t) + 1))

            self.autoloader.is_active = True
            self.currentODCalc()
        else:
            print('File not found!')

        

    def currentODCalc(self):
        self.autoloader.is_active = False

        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())

        species = IMFIT_MODES[self.mode]['Species']
        self.odK = calcOD(self.currentFile, species[0], self.mode, self.regionK)
        self.odRb = calcOD(self.currentFile, species[1], self.mode, self.regionRb)
        
        if self.fo.autoFit.isChecked() and self.frame == 'OD':
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
                
                self.currentFile.img = imageMean/float(len(x))

            self.currentODCalc()
        self.autoloader.is_active = True

    def fitCurrent(self):
        # TODO: Understand what this does and adjust to be more readable
        if self.fo.moleculeBook.isChecked():
            kAtom = 2
        else:
            kAtom = 0

        rbAtom = 1
        TOF = float(self.fo.tof.text())
        pxl = IMFIT_MODES[self.mode]["Pixel Size"]

        # TODO: Is it possible/ useful to generalize to an arbitrary number of fits with different names?
        if self.odK is not None:
            self.fitK = fitOD(self.odK, str(self.fo.kFitFunction.currentText()),kAtom,TOF,pxl)
        if self.odRb is not None:
            self.fitRb = fitOD(self.odRb, str(self.fo.rbFitFunction.currentText()), rbAtom, TOF + 6,pxl)
        self.plotCurrent()

        if self.fo.autoUpload.isChecked():
            self.process2Origin()
        print("Done fitting and uploading current shot")

    def process2Origin(self):
        imagePath = IMFIT_MODES[self.mode]["Image Path"]

        if self.fitK is not None and self.fitRb is not None:
            print("Processing K fit result")
            KProcess = processFitResult(self.fitK, imagePath)
            print("Processing Rb fit result")
            RbProcess = processFitResult(self.fitRb, imagePath)

            KProcess.data[0] = self.currentFile.fileName + "-" + imagePath
            RbProcess.data[0] = self.currentFile.fileName + "-" + imagePath

            if self.mode == 'Axial iXon Molecules In Situ': # Molecule In situ FK
                print("Uploading KRb to Origin")
                upload2Origin('KRbSpinGauss', self.fitK.fitFunction,
                                        [KProcess.data, RbProcess.data])
                return 1
            else:
                print("Uploading Rb to Origin")
                upload2Origin('Rb', self.fitRb.fitFunction, RbProcess.data)

                if self.fitK.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):
                    if self.fo.moleculeBook.isChecked():
                        KProcess.data[1] = 'KRb'
                    else:
                        KProcess.data[1] = 'K'

                    upload2Origin('K', self.fitK.fitFunction, KProcess.data)
                    return 1

                if self.fo.moleculeBook.isChecked():
                    print("Uploading KRb to Origin")
                    upload2Origin('KRb', self.fitK.fitFunction, KProcess.data)
                    return 1
                else:
                    print("Uploading K to Origin")
                    upload2Origin('K', self.fitK.fitFunction, KProcess.data)
                    return 1
            
    def plotCurrent(self):
        if self.currentFile is None:
            return
        frames = self.currentFile.frames
        species = IMFIT_MODES[self.mode]['Species']
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
            except:
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

            if self.frame == 'OD':
                image = self.odK.ODCorrected
            else:
                print(x)
                print(y)
                image = frames[species[0]][self.frame][y[0]:y[-1], x[0]:x[-1]]
            self.figs.plotUpdate(x,y,image,ch0,ch1)
    
            if self.frame == 'OD':
                if self.fitK is not None:
                    if self.fitK.fitFunction==FIT_FUNCTIONS.index('Fermi-Dirac'):
                        self.figs.plotSliceUpdate(x,[Sx,Fx],np.arange(len(R)),[R,RG,RF])
                    else:
                        self.figs.plotSliceUpdate(x,[Sx,Fx],y,[Sy,Fy])

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
            except:
                ch0 = None
                ch1 = None
                Sx = None 
                Sy = None
                Fx = None
                Fy = None
                self.figs.ax1.cla()
                self.figs.ax2.cla()

            if self.frame == 'OD':
                image = self.odRb.ODCorrected
            else:
                image = frames[species[1]][self.frame][y[0]:y[-1], x[0]:x[-1]]
            self.figs.plotUpdate(x,y,image,ch0,ch1)

            if self.frame == 'OD':
                self.figs.plotSliceUpdate(x,[Sx,Fx],y,[Sy,Fy])

    def passCamToROI(self):
        self.roi.setDefaultRegion(self.mode)

    def initializeGui(self):

        self.fo = fitOptionsWidget()
        self.figs = ImageWindows()
        self.roi = regionWidget()
        self.pf = pathWidget(self.fo, self.figs, self.roi)
        self.av = averageWidget()

        # self.pf.cameraGroup.buttonClicked.connect(self.passCamToROI)
        
        gb1 = QtWidgets.QGroupBox('File Path')
        gb1.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb1l = QtWidgets.QVBoxLayout()
        gb1l.addWidget(self.pf)
        gb1.setLayout(gb1l)

        v0 = QtWidgets.QVBoxLayout()
        v0.addStretch(4)
        v0.addWidget(gb1)
        v0.addStretch(1)

        gb2 = QtWidgets.QGroupBox('Region Selection')
        gb2.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb2l = QtWidgets.QVBoxLayout()
        h0 = QtWidgets.QHBoxLayout()
        h0.addStretch(1)
        h0.addWidget(self.roi)
        h0.addStretch(1)
        gb2l.addLayout(h0)
        gb2.setLayout(gb2l)

        v0.addWidget(gb2)
        v0.addStretch(1)
        
        gb4 = QtWidgets.QGroupBox('Averaging')
        gb4.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb4l = QtWidgets.QVBoxLayout()
        gb4l.addWidget(self.av)
        gb4.setLayout(gb4l)

        v0.addWidget(gb4)
        v0.addStretch(1)


        gb3 = QtWidgets.QGroupBox('Fitting Options')
        gb3.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
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
        self.initializeGui()

        self.autoloader = autoloader(self.pf)
        self.autoloader.connect(self.autoloader, QtCore.SIGNAL('fileArrived'), self.loadFile)

        self.autoloader.is_active = True
        self.autoloader.start()

    def saveMainImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtWidgets.QFileDialog.DontUseNativeDialog)

        ok2write = False
        if self.figs.plotTools.kSelect.isChecked():
            if hasattr(self, 'odK'):
                x = self.odK.xRange0
                y = self.odK.xRange1
                od = self.odK.ODCorrected
                ok2write = True
        elif self.figs.plotTools.rbSelect.isChecked():
            if hasattr(self, 'odRb'):
                x = self.odRb.xRange0
                y = self.odRb.xRange1
                od = self.odRb.ODCorrected
                ok2write = True
                
        if ok2write:
            f = open(xp, 'w')
            for i in range(len(y)):
                for j in range(len(x)):
                    f.write("{0:.3f},".format(od[j,i]))
                f.write("\n")
            f.close()
        else:
            print("No file to save!")

    def saveOSliceImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtWidgets.QFileDialog.DontUseNativeDialog)

        ok2write = False
        if self.figs.plotTools.kSelect.isChecked():
            if hasattr(self, 'odK') and hasattr(self.fitK.slices, 'points0'):
                x = self.odK.xRange0
                y = self.odK.xRange1
                Sx = self.fitK.slices.points0
                Fx = self.fitK.slices.fit0
                Sy = self.fitK.slices.points1
                Fy = self.fitK.slices.fit1
                ok2write = True
        elif self.figs.plotTools.rbSelect.isChecked():
            if hasattr(self, 'odRb') and hasattr(self.fitRb.slices, 'points0'):
                x = self.odRb.xRange0
                y = self.odRb.xRange1
                Sx = self.fitRb.slices.points0
                Fx = self.fitRb.slices.fit0
                Sy = self.fitRb.slices.points1
                Fy = self.fitRb.slices.fit1
                ok2write = True

        if ok2write:
            f = open(xp,'w')

            xbigger = len(x) >= len(y)
            n = max([len(x), len(y)])
            m = min([len(x), len(y)])

            for k in range(n):                
                if k < m: 
                    f.write("{0:d},{1:.3f},{2:.3f},{3:d},{4:.3f},{5:.3f}\n".format(x[k],Sx[k],Fx[k],y[k],Sy[k],Fy[k]))
                elif k > m:
                    if xbigger:
                        f.write("{0:d},{1:.3f},{2:.3f},,,\n".format(x[k],Sx[k],Fx[k]))
                    elif not xbigger:
                        f.write(",,,{0:d},{1:.3f},{2:.3f}\n".format(y[k],Sy[k],Fy[k]))
                
            f.close() 
        else:
            print("No file to save!")


    def saveRSliceImage(self):
        x = QtWidgets.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtWidgets.QFileDialog.DontUseNativeDialog)

        ok2write = False
        if self.figs.plotTools.kSelect.isChecked():
            if hasattr(self, 'odK') and hasattr(self.fitK.slices, 'radSlice'):
                R = self.fitK.slices.radSlice
                Rf = self.fitK.slices.radSliceFit
                Rg = self.fitK.slices.radSliceFitGauss
                ok2write = True
        elif self.figs.plotTools.rbSelect.isChecked():
            if hasattr(self, 'odRb') and hasattr(self.fitRb.slices, 'radSlice'):
                R = self.fitRb.slices.radSlice
                Rf = self.fitRb.slices.radSliceFit
                Rg = self.fitRb.slices.radSliceFitGauss
                ok2write = True

        if ok2write:

            f = open(xp,'w')

            for k in range(len(R)):                
                f.write("{0:.3f},{1:.3f},{2:.3f},{3:.3f}\n".format(k,R[k],Rf[k],Rg[k]))
                
            f.close() 
        else:
            print("No file to save!")

    def loadFromMenu(self): 

        d = str(self.pf.filePath.text())
        if self.pf.cameraGroup.checkedId():
            ext = '*.csv'
        else:
            ext = '*.dat'

        x = QtWidgets.QFileDialog()
        xp = x.getOpenFileName(self,'Select a file to load', filter=ext, directory=d,options=QtWidgets.QFileDialog.DontUseNativeDialog)

        self.pf.filePath.setText(xp)
        self.loadFile()

    def getStyleSheet(self,path):
        f = QtCore.QFile(path)
        f.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text)
        stylesheet = QtCore.QTextStream(f).readAll()
        f.close()

        return stylesheet

    def createToolbar(self):

        exitAction = QtWidgets.QAction("Exit",self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(QtWidgets.qApp.quit)

        loadAction = QtWidgets.QAction("Load Image",self)
        loadAction.setShortcut('Ctrl+O')
        loadAction.triggered.connect(self.loadFromMenu)

        refreshAction = QtWidgets.QAction("Refresh all", self)
        refreshAction.setShortcut('Ctrl+R')
        refreshAction.triggered.connect(self.refreshGui)

        saveMain = QtWidgets.QAction("Save Main Image", self)
        saveMain.setShortcut('Ctrl+S')
        saveMain.triggered.connect(self.saveMainImage)
        saveOSlice = QtWidgets.QAction("Save Ortho Slice", self)
        saveOSlice.triggered.connect(self.saveOSliceImage)
        saveRSlice = QtWidgets.QAction("Save Radial Slice", self)
        saveRSlice.triggered.connect(self.saveRSliceImage)


        menubar = self.menuBar()
        
        fileMenu = menubar.addMenu("File")
        fileMenu.addAction(loadAction)
        fileMenu.addAction(refreshAction)
        fileMenu.addAction(exitAction)

        saveMenu = menubar.addMenu("Save")
        saveMenu.addAction(saveMain)
        saveMenu.addAction(saveOSlice)
        saveMenu.addAction(saveRSlice)
        

    def closeEvent(self, event):
        self.autoloader.terminate()

if __name__=="__main__":
    # The following two lines tell windows that python is only hosting this application
    myappid = 'krb.imfitdue' # arbitrary string
    if os.name == 'nt':
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    
    app = QtWidgets.QApplication(sys.argv)
    w = imfitDue()
    w.setGeometry(100,100, 1200, 800)

    appico = QtGui.QIcon('main.ico')
    w.setWindowIcon(appico)

    w.show()
    sys.exit(app.exec_())

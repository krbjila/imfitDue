import sys, os
sys.path.append('./lib')

from PyQt4 import QtGui, QtCore

from imageRead import *
from imageProcess import *
from imfitDefaults import *
from gui_helpers import *

import numpy as np

class imfitDue(QtGui.QMainWindow):
    def __init__(self,Parent=None):
        super(imfitDue, self).__init__(Parent)

        self.regionRb = [0]*4
        self.regionK = [0]*4


        self.initializeGui()
        self.createToolbar()
        self.makeConnections()

        self.autoloader = autoloader(self.pf)
        self.autoloader.connect(self.autoloader, QtCore.SIGNAL('fileArrived'), self.loadFile)
        self.autoloader.start()

        self.fitK = None
        self.fitRb = None

    def makeConnections(self):
        self.pf.loadButton.clicked.connect(self.loadFile)
        self.fo.fitButton.clicked.connect(self.fitCurrent)

        self.figs.plotTools.atomSelectGroup.buttonClicked.connect(self.plotCurrent)

        for i in range(2):
            for j in range(4):
                self.roi.region[i][j].returnPressed.connect(self.currentODCalc)
        
        self.av.averageButton.clicked.connect(self.averageImages)
        self.fo.uploadButton.clicked.connect(self.process2Origin)

    def loadFile(self):
        path = str(self.pf.filePath.text())
        if os.path.isfile(path):
            # Terminate the autoloader thread so that the reader doesn't get bogged down
            self.autoloader.terminate()
            if self.pf.cameraGroup.checkedId() == 0:
                # self.currentFile = readPIImage(path) # Changed to xiQ on 9/26/2019
                self.currentFile = readXimeaImage(path)
            elif self.pf.cameraGroup.checkedId() == 1:
                self.currentFile = readIXONImage(path)
            

            if self.pf.autoLoad.isChecked():
                t = self.pf.autoLoadFile.text()
                self.pf.autoLoadFile.setText(str(int(t) + 1))

            self.autoloader.start()

        else:
            print('File not found!')

        self.currentODCalc()

    def currentODCalc(self):
        self.autoloader.terminate()
        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())


        imagePath = IMAGING_PATHS.index(self.fo.imagePath.currentText())
        self.odK = calcOD(self.currentFile,'K', imagePath,self.regionK)
        self.odRb = calcOD(self.currentFile, 'Rb', imagePath, self.regionRb)

        
        if self.fo.autoFit.isChecked():
            self.fitCurrent()
        else:
            self.fitK = None
            self.fitRb = None
            self.plotCurrent()
        self.autoloader.start()

    def averageImages(self):
        self.autoloader.terminate()
        self.pf.autoLoad.setChecked(False)
        x = self.av.getFileNumbers()
        p = self.pf.cameraGroup.checkedId() 

        if p == 0:
            # path = DEFAULT_PATH_PI + "pi_{}.spe"  # Changed to xiQ on 9/26/2019
            path = DEFAULT_PATH_XIMEA + "xi_{}.dat"
        elif p == 1:
            path = DEFAULT_PATH_IXON + "ixon_{}.csv"

        firstFile = True
        if x is not None:
            for k in x:
                print(path.format(k))
                if p == 0:
                    # self.currentFile = readPIImage(path.format(k)) # Changed to xiQ on 9/26/2019
                    self.currentFile = readXimeaImage(path.format(k))
                    if firstFile:
                        imageMean = self.currentFile.img
                        firstFile = False
                    else:
                        imageMean += self.currentFile.img
                    
                    self.currentFile.img = imageMean/float(len(x))
                elif p == 1:
                    self.currentFile = readIXONImage(path.format(k))
                    if firstFile:
                        imageMean = self.currentFile.img
                        firstFile = False
                    else:
                        imageMean += self.currentFile.img
                    
                    self.currentFile.img = imageMean/float(len(x))

            self.currentODCalc()
            self.autoloader.start()



    def fitCurrent(self):
       

        if self.fo.moleculeBook.isChecked():
            kAtom = 2
        else:
            kAtom = 0

        rbAtom = 1
        TOF = float(self.fo.tof.text())
        pxl = PIXEL_SIZES[IMAGING_PATHS.index(self.fo.imagePath.currentText())]

        self.fitK = fitOD(self.odK, str(self.fo.kFitFunction.currentText()),kAtom,TOF,pxl)
        self.fitRb = fitOD(self.odRb, str(self.fo.rbFitFunction.currentText()), rbAtom, TOF + 6,pxl)
        self.plotCurrent()

        if self.fo.autoUpload.isChecked():
            self.process2Origin()

    def process2Origin(self):
        
        imagePath = self.fo.imagePath.currentText()

        if self.fitK is not None and self.fitRb is not None:
            KProcess = processFitResult(self.fitK, imagePath)
            RbProcess = processFitResult(self.fitRb, imagePath)

            KProcess.data[0] = self.currentFile.fileName[0] + "-" + IMAGING_PATHS[KProcess.imagePath]
            RbProcess.data[0] = self.currentFile.fileName[0] + "-" + IMAGING_PATHS[KProcess.imagePath]



            upload2Origin('Rb', self.fitRb.fitFunction, RbProcess.data)

            if self.fitK.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):
                if self.fo.moleculeBook.isChecked():
                    KProcess.data[1] = 'KRb'
                else:
                    KProcess.data[1] = 'K'

                upload2Origin('K', self.fitK.fitFunction, KProcess.data)
                return 1


            if self.fo.moleculeBook.isChecked():
                upload2Origin('KRb', self.fitK.fitFunction, KProcess.data)
                return 1
            else:
                upload2Origin('K', self.fitK.fitFunction, KProcess.data)
                return 1
            






    def plotCurrent(self):
        
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

            self.figs.plotUpdate(x,y,self.odK.ODCorrected,ch0,ch1)

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

            self.figs.plotUpdate(x,y,self.odRb.ODCorrected,ch0,ch1)
            self.figs.plotSliceUpdate(x,[Sx,Fx],y,[Sy,Fy])

    def passCamToROI(self):
        # iXon is 1, XIMEA is 0
        if self.pf.cameraGroup.checkedId():
            self.roi.setDefaultRegion('IXON')
        else:
            self.roi.setDefaultRegion('XIMEA')


    def initializeGui(self):

        self.pf = pathWidget()
        self.roi = regionWidget()
        self.fo = fitOptionsWidget()
        self.av = averageWidget()

        self.pf.cameraGroup.buttonClicked.connect(self.passCamToROI)


        self.figs = ImageWindows()
        
        gb1 = QtGui.QGroupBox('File Path')
        gb1.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb1l = QtGui.QVBoxLayout()
        gb1l.addWidget(self.pf)
        gb1.setLayout(gb1l)

        v0 = QtGui.QVBoxLayout()
        v0.addStretch(4)
        v0.addWidget(gb1)
        v0.addStretch(1)

        gb2 = QtGui.QGroupBox('Region Selection')
        gb2.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb2l = QtGui.QVBoxLayout()
        h0 = QtGui.QHBoxLayout()
        h0.addStretch(1)
        h0.addWidget(self.roi)
        h0.addStretch(1)
        gb2l.addLayout(h0)
        gb2.setLayout(gb2l)

        v0.addWidget(gb2)
        v0.addStretch(1)
        
        gb4 = QtGui.QGroupBox('Averaging')
        gb4.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb4l = QtGui.QVBoxLayout()
        gb4l.addWidget(self.av)
        gb4.setLayout(gb4l)

        v0.addWidget(gb4)
        v0.addStretch(1)


        gb3 = QtGui.QGroupBox('Fitting Options')
        gb3.setStyleSheet(self.getStyleSheet('./lib/styles.qss'))
        gb3l = QtGui.QVBoxLayout()
        gb3l.addWidget(self.fo)
        gb3.setLayout(gb3l)

        v0.addWidget(gb3)

        h = QtGui.QHBoxLayout()
        h.addLayout(v0)
        h.addWidget(self.figs)

        self.mainWidget = QtGui.QWidget()
        self.mainWidget.setAutoFillBackground(True)
        p = self.mainWidget.palette()
        p.setColor(self.mainWidget.backgroundRole(), QtCore.Qt.white)
        self.mainWidget.setPalette(p)
        self.mainWidget.setLayout(h)

        self.setCentralWidget(self.mainWidget)

    def refreshGui(self):
        self.autoloader.terminate()
        self.initializeGui()

        self.autoloader = autoloader(self.pf)
        self.autoloader.connect(self.autoloader, QtCore.SIGNAL('fileArrived'), self.loadFile)
        self.autoloader.start()

    def saveMainImage(self):
        x = QtGui.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtGui.QFileDialog.DontUseNativeDialog)

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
        x = QtGui.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtGui.QFileDialog.DontUseNativeDialog)

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
        x = QtGui.QFileDialog()
        xp = x.getSaveFileName(self, "Save image as", "untitled.dat", "Data file (*.dat)",options=QtGui.QFileDialog.DontUseNativeDialog)

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
            # ext = '*.spe'
            ext = '*.dat'

        x = QtGui.QFileDialog()
        xp = x.getOpenFileName(self,'Select a file to load', filter=ext, directory=d,options=QtGui.QFileDialog.DontUseNativeDialog)

        self.pf.filePath.setText(xp)
        self.loadFile()


    def getStyleSheet(self,path):
        f = QtCore.QFile(path)
        f.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text)
        stylesheet = QtCore.QTextStream(f).readAll()
        f.close()

        return stylesheet

    def createToolbar(self):

        exitAction = QtGui.QAction("Exit",self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(QtGui.qApp.quit)

        loadAction = QtGui.QAction("Load Image",self)
        loadAction.setShortcut('Ctrl+O')
        loadAction.triggered.connect(self.loadFromMenu)

        refreshAction = QtGui.QAction("Refresh all", self)
        refreshAction.setShortcut('Ctrl+R')
        refreshAction.triggered.connect(self.refreshGui)

        saveMain = QtGui.QAction("Save Main Image", self)
        saveMain.setShortcut('Ctrl+S')
        saveMain.triggered.connect(self.saveMainImage)
        saveOSlice = QtGui.QAction("Save Ortho Slice", self)
        saveOSlice.triggered.connect(self.saveOSliceImage)
        saveRSlice = QtGui.QAction("Save Radial Slice", self)
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

    app = QtGui.QApplication(sys.argv)

    w = imfitDue()
    w.setGeometry(100,100, 1200, 800)
    w.show()

    sys.exit(app.exec_())

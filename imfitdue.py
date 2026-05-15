# TODO: use image metadata to automatically set binning and change region selection
# TODO: add support for an arbitrary number of images, not just K and Rb frames

import sys, os

from PyQt5 import QtWidgets, QtCore, QtGui
import ctypes

from lib.imageRead import *
from lib.imageProcess import *
from lib.imfitDefaults import *
from lib.gui_helpers import *

import numpy as np

from pymongo import MongoClient, errors
from bson.json_util import loads, dumps

from datetime import datetime

# ADDED FOR DEBUGGING PURPOSES
class PlotWindow(QtWidgets.QMainWindow):
    def __init__(self, img):
        super().__init__()
        self.setWindowTitle("Matplotlib Window")

        fig = Figure()
        canvas = FigureCanvas(fig)
        self.setCentralWidget(canvas)

        ax = fig.add_subplot(111)
        ax.imshow(img, cmap='viridis')
        ax.set_title("My Plot")

        canvas.draw()

def show_warning_messagebox_defringing(warning_text="Warning: Have you initialized the defringing? \nPerhaps the region changed?"):
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Icon.Warning)
    msg.setText(warning_text)
    msg.setWindowTitle("Warning MessageBox")
    msg.setStandardButtons(QtWidgets.QMessageBox.StandardButton.Ok | QtWidgets.QMessageBox.StandardButton.Cancel)
    retval = msg.exec()

class imfitDue(QtWidgets.QMainWindow):
    def __init__(self, Parent=None):
        super(imfitDue, self).__init__(Parent)

        self.setWindowTitle("ImfitDue: KRb Image Fitting")

        self.regionRb = [0] * 4
        self.regionK = [0] * 4
        self.pRb = [0] * 4
        self.pK = [0] * 4
        self.fK = [0] * 2

        self.sigblur = 0

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
                
        self.Ref_Mat_Full_lgK = None
        self.Ref_Mat_lgK = None
        self.Ref_Mat_Full_shK = None
        self.Ref_Mat_shK = None
        self.RFRFT_inv_lgK = None
        self.RFRFT_inv_shK = None
        self.ind_maskK = None

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
                self.roi.region[i][j].returnPressed.connect(self.updateCrop)
                
        for i in range(2):
            for j in range(2):
                self.roi.fitregion[i][j].returnPressed.connect(self.cropaverageImages)

        self.av.averageButton.clicked.connect(self.averageImages)
        self.av.initdefrButton.clicked.connect(self.intializeDefringe)
        self.av.gaussblurButton.clicked.connect(self.GaussBlur)
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

    def currentBGCalc(self):
        self.autoloader.is_active = False

        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())

        species = IMFIT_MODES[self.mode]["Species"]
        try:
            self.odKBG = calcOD(self.BGFile, species[0], self.mode, self.regionK, self.sigblur)
            self.odRbBG = calcOD(self.BGFile, species[1], self.mode, self.regionRb, self.sigblur)
        except Exception as e:
            print("Could not calculate BG: {}".format(e))

        self.autoloader.is_active = True
        print("Done calculating current BG")

    def currentODCalc(self):
        self.autoloader.is_active = False

        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())
            self.pK[i] = float(self.av.region_p[0][i].text())
            self.pRb[i] = float(self.av.region_p[1][i].text())
            
        for i in range(2):
            self.fK[i] = float(self.roi.fitregion[0][i].text())

        species = IMFIT_MODES[self.mode]["Species"]
        try:
            if self.av.b1_nobg.isChecked():
                self.odK = calcOD(self.currentFile, species[0], self.mode, self.regionK, self.sigblur)
                self.odRb = calcOD(self.currentFile, species[1], self.mode, self.regionRb, self.sigblur)

            elif self.av.b2_bgsu.isChecked():
                self.odK = calcOD(self.currentFile, species[0], self.mode, self.regionK, self.sigblur)
                self.odRb = calcOD(self.currentFile, species[1], self.mode, self.regionRb, self.sigblur)
                self.odK.OD = self.odK.OD - self.odKBG.OD
                self.odRb.OD = self.odRb.OD - self.odRbBG.OD
                self.odK.ODCorrected = self.odK.ODCorrected - self.odKBG.ODCorrected
                self.odRb.ODCorrected = self.odRb.ODCorrected - self.odRbBG.ODCorrected
                self.odK.n = self.odK.n - self.odKBG.n
                self.odRb.n = self.odRb.n - self.odRbBG.n
                self.odK.nerr = np.sqrt(self.odK.nerr**2 + self.odKBG.nerr**2)
                self.odRb.nerr = np.sqrt(self.odRb.nerr**2 + self.odRbBG.nerr**2)

            elif self.av.b3_defr.isChecked():
                # The plan is to defringe in the averaging function, so here we just calculate the OD normally
                self.odK = calcOD(self.currentFile, species[0], self.mode, self.regionK, self.sigblur)
                self.odRb = calcOD(self.currentFile, species[1], self.mode, self.regionRb, self.sigblur)

                # fitbox = self.fK
                fitregionK = self.regionK.copy()
                fitregionK[2] = self.fK[0]
                fitregionK[3] = self.fK[1]
                self.odKfit = calcOD(self.currentFile, species[0], self.mode, fitregionK, self.sigblur)
                # For debugging purposes only:
                # self.plot_window = PlotWindow(self.odKfit.ODCorrected)
                # self.plot_window.show()

        except Exception as e:
            print("Could not calculate OD: {}".format(e))

        if (
            self.fo.autoFit.isChecked()
            and self.frame == "OD"
            or self.frame == "Column Density"
        ):
            self.fitCurrent()
            print("Done fitting current image")
        else:
            self.fitK = None
            self.fitRb = None
            self.plotCurrent()

        self.autoloader.is_active = True
        print("Done calculating current OD")

    def updateCrop(self):
        for i in range(2):
            for j in range(2):
                if float(self.roi.fitregion[i][j].text()) > float(self.roi.region[i][j + 2].text()):
                    self.roi.fitregion[i][j].setText(self.roi.region[i][j + 2].text())

        if self.av.b2_bgsu.isChecked():
            self.averageImages()
        else:
            self.currentODCalc()
        # also set the defringe particle region xc and yc to the crop region xc and yc
        for i in range(2):
           self.av.region_p[0][i].setText(self.roi.region[0][i].text())
           self.av.region_p[1][i].setText(self.roi.region[1][i].text())


    def intializeDefringe(self):  # FOR NOW IMPLEMENTING ONLY FOR iXon side and for K and Rb
        # At the moment this code is quite messy: I copy pasted a lot of code from elsewhere. It could have been a function...
        for i in range(4):
            self.regionK[i] = float(self.roi.region[0][i].text())
            self.regionRb[i] = float(self.roi.region[1][i].text())
            self.pK[i] = float(self.av.region_p[0][i].text())
            self.pRb[i] = float(self.av.region_p[1][i].text())

        particle_regK = [self.pK[0], self.pK[1], self.pK[2], self.pK[3]]
        particle_regRb = [self.pRb[0], self.pRb[1], self.pRb[2], self.pRb[3]]

        regionK = [self.regionK[0], self.regionK[1], self.regionK[2], self.regionK[3]]
        regionRb = [self.regionRb[0], self.regionRb[1], self.regionRb[2], self.regionRb[3]]

        x_recK = int(particle_regK[0] - particle_regK[2]/2 - (regionK[0] - regionK[2]/2))
        y_recK = int(particle_regK[1] - particle_regK[3]/2 - (regionK[1] - regionK[3]/2))

        x_recRb = int(particle_regRb[0] - particle_regRb[2]/2 - (regionRb[0] - regionRb[2]/2))
        y_recRb = int(particle_regRb[1] - particle_regRb[3]/2 - (regionRb[1] - regionRb[3]/2))

        maskK = np.zeros((int(regionK[3]), int(regionK[2]))) + 1
        maskK[
            y_recK : y_recK + int(particle_regK[3]),
            x_recK : x_recK + int(particle_regK[2]),
        ] = 0
        ind_maskK = np.where(maskK.ravel() != 0)
        self.ind_maskK = ind_maskK

        y = self.av.getBackgroundFileNumbers()

        path = str(
                self.av.bgPath.text()
            )

        species = IMFIT_MODES[self.mode]["Species"]

        sz_frame = 0
        # Use first image to get size of full light frames
        if y is not None:
            try:
                initializer = readImage(
                    self.mode, (path + IMFIT_MODES[self.mode]["Default Suffix"]).format(y[0])
                )
                # # DEBUG
                # test = initializer.getFrame(species[0], "Shadow")
                # self.plot_window = PlotWindow(test[
                #                                 y_recK : y_recK + int(particle_regK[3]),
                #                                 x_recK : x_recK + int(particle_regK[2]),
                #                             ])
                # self.plot_window.show()

                if initializer is None:
                    return
                
                sz_frame = np.size(initializer.getFrame(species[0], "Shadow"))

            except Exception as e:
                print("Could not defringe images: {}".format(e))
        
        Ref_Mat_Full_lgK = np.zeros((sz_frame, np.size(y)))
        Ref_Mat_Full_shK = np.zeros((sz_frame, np.size(y)))

        Ref_Mat_lgK = np.zeros((np.size(ind_maskK), np.size(y)))
        Ref_Mat_shK = np.zeros((np.size(ind_maskK), np.size(y)))


        r0K = int(regionK[0] - np.floor(regionK[2] / 2))
        r1K = int(regionK[0] + np.floor(regionK[2] / 2))
        r2K = int(regionK[1] - np.floor(regionK[3] / 2))
        r3K = int(regionK[1] + np.floor(regionK[3] / 2))

        xRange0K = range(r0K, r1K)
        xRange1K = range(r2K, r3K)
        
        if y is not None:
            try:
                i_cnt = 0
                for k in y:
                    print(
                        "Reading file: {}".format(
                            (path + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)
                        )
                    )
                    self.BGFile = readImage(
                        self.mode, (path + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)
                    )
                    if self.BGFile is None:
                        return
                    shadow = self.BGFile.getFrame(species[0], "Shadow")
                    light = self.BGFile.getFrame(species[0], "Light")
                    dark = self.BGFile.getFrame(species[0], "Dark")

                    shadowCrop = cropArray(shadow, xRange1K, xRange0K)
                    lightCrop = cropArray(light, xRange1K, xRange0K)
                    darkCrop = cropArray(dark, xRange1K, xRange0K)

                    s1 = shadowCrop - darkCrop
                    s2 = lightCrop - darkCrop
                    s1f = shadow - dark
                    s2f = light - dark
                                    
                    Ref_Mat_Full_lgK[:, i_cnt] = s2f.ravel()
                    Ref_Mat_lgK[:, i_cnt] = s2.ravel()[ind_maskK]

                    Ref_Mat_Full_shK[:, i_cnt] = s1f.ravel()
                    Ref_Mat_shK[:, i_cnt] = s1.ravel()[ind_maskK]
                    i_cnt += 1

                # TO DO: Return these matrices to the main imfitDue class
                # then use them in computing the corrected OD
                self.Ref_Mat_Full_lgK = Ref_Mat_Full_lgK
                self.Ref_Mat_Full_shK = Ref_Mat_Full_shK

                Ref_Mat_lgK = Ref_Mat_lgK.T
                RFRFT_inv_lgK = np.linalg.inv(Ref_Mat_lgK @ Ref_Mat_lgK.T)
                Ref_Mat_shK = Ref_Mat_shK.T
                RFRFT_inv_shK = np.linalg.inv(Ref_Mat_shK @ Ref_Mat_shK.T)

                self.Ref_Mat_lgK = Ref_Mat_lgK
                self.RFRFT_inv_lgK = RFRFT_inv_lgK
                self.Ref_Mat_shK = Ref_Mat_shK
                self.RFRFT_inv_shK = RFRFT_inv_shK

            except Exception as e:
                print("Could not defringe images: {}".format(e))
    
    def GaussBlur(self):
        sigblur_str = self.av.sigblur.text()
        try:
            sigblur = float(sigblur_str)
            if sigblur < 0:
                raise ValueError("Signal blur value must be zero or positive.")
            else:
                self.sigblur = sigblur
                self.averageImages()
        except ValueError as e:
            print("Invalid signal blur value: {}".format(e))
            show_warning_messagebox_defringing("Invalid signal blur value: {}".format(e))
            return

    def cropaverageImages(self):
        update_flag = True
        # Check if the selected fitting range lies within the overall crop region
        for i in range(2):
            for j in range(2):
                if float(self.roi.fitregion[i][j].text()) > float(self.roi.region[i][j + 2].text()):
                    update_flag = False
        if update_flag:
            # self.intializeDefringe()
            self.averageImages()
        else:
            warn_str = 'Please select a fitting range that is smaller than the cropped region!'
            show_warning_messagebox_defringing(warn_str)

    def averageImages(self):  # FOR NOW IMPLEMENTING ONLY FOR iXon Side
        # For now the background subtraction has not been extensively tested yet
        self.autoloader.is_active = False

        self.pf.autoLoad.setChecked(False)
        try:
            x = self.av.getFileNumbers()
            print(str(self.pf.filePath.text()))

            path = str(
                self.pf.filePath.text()
            )  # IMFIT_MODES[self.mode]["Default Path"]
            bgpath = str(
                self.av.bgPath.text()
            ) 

            defringe_flag = False
            if self.av.b3_defr.isChecked():
                if np.any(self.Ref_Mat_Full_lgK == None):
                    show_warning_messagebox_defringing()
                else:
                    print("Defringing during averaging is not yet implemented for Rb.")
                    defringe_flag = True

            if self.av.b2_bgsu.isChecked():
                y = self.av.getBackgroundFileNumbers()
                # self.BGFile = None # reinitialize BGFiles - didn't resolve issue
                # self.odKBG = None
                # self.odRbBG = None
                firstBGFile = True
                if y is not None:
                    for k in y:
                        
                        print(
                            "Subtracting file: {}".format(
                                (bgpath + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)
                            )
                        )
                        self.BGFile = readImage(
                            self.mode, (bgpath + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)
                        )
                        if self.BGFile is None:
                            return
                        if firstBGFile:
                            bg_frame_dict = self.BGFile.frames
                            species_list = list(self.BGFile.frames.keys())
                            firstBGFile = False
                        else:
                            species_list = list(self.BGFile.frames.keys())
                            for idx, species in enumerate(species_list):
                                frame_list = list(self.BGFile.frames[species].keys())
                                for idy, frame_name in enumerate(frame_list):
                                    species_frame_bg = self.BGFile.frames[species][
                                        frame_name
                                    ]
                                    bg_frame_dict[species][
                                        frame_name
                                    ] += species_frame_bg

                    for idx, species in enumerate(species_list):
                        frame_list = list(self.BGFile.frames[species].keys())
                        for idy, frame_name in enumerate(frame_list):
                            # Actually put the average image in the current image dict
                            self.BGFile.frames[species][frame_name] = bg_frame_dict[
                                species
                            ][frame_name] / float(len(y))
                    self.currentBGCalc()

            firstFile = True
            if x is not None:
                for k in x:
                    print("Loading file: {}".format((path + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)))
                    self.currentFile = readImage(
                        self.mode, (path + IMFIT_MODES[self.mode]["Default Suffix"]).format(k)
                    )
                    if self.currentFile is None:
                        return
                    if firstFile:
                        avg_frame_dict = self.currentFile.frames
                        species_list = list(self.currentFile.frames.keys())
                        
                        if defringe_flag:
                            regionK = [self.regionK[0], self.regionK[1], self.regionK[2], self.regionK[3]]
                            r0K = int(regionK[0] - np.floor(regionK[2] / 2))
                            r1K = int(regionK[0] + np.floor(regionK[2] / 2))
                            r2K = int(regionK[1] - np.floor(regionK[3] / 2))
                            r3K = int(regionK[1] + np.floor(regionK[3] / 2))

                            xRange0K = range(r0K, r1K)
                            xRange1K = range(r2K, r3K)

                            light0 = avg_frame_dict["K"]["Light"] - avg_frame_dict["K"]["Dark"]
                            light1 = cropArray(light0, xRange1K, xRange0K)
                            shad0 = avg_frame_dict["K"]["Shadow"] - avg_frame_dict["K"]["Dark"]
                            shad1 = cropArray(shad0, xRange1K, xRange0K)

                            reg_sz_np = self.regionK[2]*self.regionK[3] - self.pK[2]*self.pK[3]

                            if reg_sz_np != np.size(self.Ref_Mat_lgK[0,:]):
                                show_warning_messagebox_defringing()


                            print(np.size(self.Ref_Mat_lgK[0,:]))
                            print(np.size(xRange0K)*np.size(xRange1K))
                            print(np.size(light1.ravel()[self.ind_maskK]))

                            a_F_lg = light1.ravel()[self.ind_maskK].T # following Vogel notation
                            a_F_sh = shad1.ravel()[self.ind_maskK].T  # following Vogel notation
                            
                            t_lg = self.Ref_Mat_lgK @ a_F_lg
                            weight_vec_lg = self.RFRFT_inv_lgK @ t_lg

                            t_sh = self.Ref_Mat_shK @ a_F_sh
                            weight_vec_sh = self.RFRFT_inv_shK @ t_sh

                            a_cor_lg = light0.ravel() - weight_vec_lg @ self.Ref_Mat_Full_lgK.T
                            light0_defringe = a_cor_lg.reshape(np.shape(light0)) + light0

                            a_cor_sh = shad0.ravel() - weight_vec_sh @ self.Ref_Mat_Full_shK.T
                            shad0_defringe = a_cor_sh.reshape(np.shape(light0)) + light0

                            # Dirty solution - since defringing was computed with dark frames already subtracted
                            # we just set the dark frame to zero in the following
                            avg_frame_dict["K"]["Light"] = light0_defringe
                            avg_frame_dict["K"]["Shadow"] = shad0_defringe
                            avg_frame_dict["K"]["Dark"] = light0_defringe*0

                        firstFile = False
                    else:
                        frames, metadata = self.currentFile.getData()
                        species_list = list(self.currentFile.frames.keys())
                        for idx, species in enumerate(species_list):

                            if defringe_flag and species == "K":
                                    light0 = self.currentFile.frames["K"]["Light"] - self.currentFile.frames["K"]["Dark"]
                                    light1 = cropArray(light0, xRange1K, xRange0K)
                                    shad0 = self.currentFile.frames["K"]["Shadow"] - self.currentFile.frames["K"]["Dark"]
                                    shad1 = cropArray(shad0, xRange1K, xRange0K)

                                    a_F_lg = light1.ravel()[self.ind_maskK].T  # following Vogel notation
                                    a_F_sh = shad1.ravel()[self.ind_maskK].T  # following Vogel notation
                                    
                                    t_lg = self.Ref_Mat_lgK @ a_F_lg
                                    weight_vec_lg = self.RFRFT_inv_lgK @ t_lg

                                    t_sh = self.Ref_Mat_shK @ a_F_sh
                                    weight_vec_sh = self.RFRFT_inv_shK @ t_sh

                                    a_cor_lg = light0.ravel() - weight_vec_lg @ self.Ref_Mat_Full_lgK.T
                                    light0_defringe = a_cor_lg.reshape(np.shape(light0)) + light0

                                    a_cor_sh = shad0.ravel() - weight_vec_sh @ self.Ref_Mat_Full_shK.T
                                    shad0_defringe = a_cor_sh.reshape(np.shape(light0)) + light0

                                    # Dirty solution - since defringing was computed with dark frames already subtracted
                                    # we just set the dark frame to zero in the following
                                    avg_frame_dict["K"]["Light"] += light0_defringe
                                    avg_frame_dict["K"]["Shadow"] += shad0_defringe
                                    avg_frame_dict["K"]["Dark"] += light0_defringe*0

                            else:
                                frame_list = list(self.currentFile.frames[species].keys())
                                for idy, frame_name in enumerate(frame_list):
                                    species_frame = self.currentFile.frames[species][
                                        frame_name
                                    ]
                                    avg_frame_dict[species][frame_name] += species_frame

                for idx, species in enumerate(species_list):
                    frame_list = list(self.currentFile.frames[species].keys())
                    for idy, frame_name in enumerate(frame_list):
                        # Actually put the average image in the current image dict and divide by the number of imgs
                        self.currentFile.frames[species][frame_name] = avg_frame_dict[
                            species
                        ][frame_name] / float(len(x))

                self.currentODCalc()
                self.currentFile.fileName = "Average of " + str(x)
                print(str(self.currentFile.fileName))
            self.autoloader.is_active = True
        except Exception as e:
            print("Could not average images: {}".format(e))

    def fitCurrent(self):
        # TODO: Understand what this does and adjust to be more readable
        rbAtom = 1
        kAtom = 0
        TOF = 0
        fx = 0
        fy = 0
        fz = 0
        Nscaler = 1 #fudge factor for the number of particles in the fiex betamu fit

        if self.fo.tof.text() != "":
            TOF = float(self.fo.tof.text())
        if self.fo.fx.text() != "":
            fx = float(self.fo.fx.text())
        if self.fo.fy.text() != "":
            fy = float(self.fo.fy.text())
        if self.fo.fz.text() != "":
            fz = float(self.fo.fz.text())
        if self.fo.Nscaler.text() != "":
            Nscaler = float(self.fo.Nscaler.text())
        
        mass_betamu_fit_species = "K"
        if self.fo.kMass.isChecked():
            mass_betamu_fit_species = "K"
        if self.fo.krbMass.isChecked():
            mass_betamu_fit_species = "KRb"

        WingRad = 0.0
        if self.fo.gausswingrad.text() != "":
            WingRad = float(self.fo.gausswingrad.text())

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
                if self.av.b3_defr.isChecked():
                    # fitbox = self.fK
                    # So here we'll be checking whether the fit box fits within the inside of the plotbox
                    # Or should we just do it for plotting and defringing?

                    self.fitK = fitOD(
                        self.mode,
                        self.odKfit,
                        str(self.fo.kFitFunction.currentText()),
                        kAtom,
                        TOF,
                        pxl,
                        WingRad = WingRad, # Added WingRad parameter to fitOD for Gaussian Wing Fitting
                        fx = fx,
                        fy = fy,
                        fz = fz,
                        mbemu = mass_betamu_fit_species,
                        Nscaler = Nscaler,
                        fullOD = self.odK,
                    )
                else:
                    self.fitK = fitOD(
                        self.mode,
                        self.odK,
                        str(self.fo.kFitFunction.currentText()),
                        kAtom,
                        TOF,
                        pxl,
                        WingRad = WingRad, # Added WingRad parameter to fitOD for Gaussian Wing Fitting
                        fx = fx,
                        fy = fy,
                        fz = fz,
                        mbemu = mass_betamu_fit_species,
                        Nscaler = Nscaler,
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
                    WingRad=WingRad,  # Added WingRad parameter to fitOD for Gaussian Wing Fitting
                    fx = fx,
                    fy = fy,
                    fz = fz,
                    mbemu = mass_betamu_fit_species,
                    Nscaler = Nscaler,
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
                if "Gauss" in FIT_FUNCTIONS[KProcess.fitObject.fitFunction]:
                    if "Mask" in FIT_FUNCTIONS[KProcess.fitObject.fitFunction]:
                        print("Uploading KRb to Origin")
                        upload2Origin(
                            "KRbSpinGaussMask",
                            self.fitK.fitFunction,
                            [KProcess.data, RbProcess.data],
                        )
                        print("Done uploading KRb to Origin")
                    else:
                        print("Uploading KRb to Origin")
                        upload2Origin(
                            "KRbSpinGauss",
                            self.fitK.fitFunction,
                            [KProcess.data, RbProcess.data],
                        )
                        print("Done uploading KRb to Origin")
                else:
                    print("Uploading KRb to Origin")
                    upload2Origin("N0", self.fitK.fitFunction, KProcess.data)
                    upload2Origin("N1", self.fitRb.fitFunction, RbProcess.data)
                    print("Done uploading Fermi-Dirac KRb to Origin")
                return 1
            else:
                if "Gauss" not in FIT_FUNCTIONS[KProcess.fitObject.fitFunction]:
                    print("Uploading KRb to Origin")
                    upload2Origin("KRb", self.fitK.fitFunction, KProcess.data)
                    print("Done uploading Fermi-Dirac KRb to Origin")
                else:
                    print("Uploading KRb to Origin")
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

            box = self.pK if self.av.b3_defr.isChecked() else None
            fitbox = self.fK if self.av.b3_defr.isChecked() else None
            print("Fit box:")
            print(fitbox)
            
            x = self.odK.xRange0
            y = self.odK.xRange1

            try:
                ch0 = self.fitK.slices.ch0
                ch1 = self.fitK.slices.ch1
                Sx = self.fitK.slices.points0
                Sy = self.fitK.slices.points1
                Fx = self.fitK.slices.fit0
                Fxa = self.fitK.slices.fit0a
                Fy = self.fitK.slices.fit1
                Fya = self.fitK.slices.fit1a

                R = self.fitK.slices.radSlice
                RG = self.fitK.slices.radSliceFitGauss
                RF = self.fitK.slices.radSliceFit
                if box is None: # old box usage below, not sure what the old usage of the box was
                    if hasattr(self.fitK, "box"):
                        box = self.fitK.box
                    else:
                        box = None
            
            except Exception as e:
                ch0 = None
                ch1 = None
                Sx = None
                Sy = None
                Fx = None
                Fxa = None
                Fy = None
                Fya = None
                # box = None

                # Initialized as 'None' in case the fit function does not have radial slices
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
            self.figs.plotUpdate(x, y, image, ch0, ch1, box, fitbox)

            # Adjust slices to fitbox if defined
            if fitbox is not None:
                nx = fitbox[0]//2*2  # ensure even number (ch0 and ch1 always even length)
                start = (len(x) - nx) // 2
                x = x[int(start):int(start + nx)]

                ny = fitbox[1]//2*2  # ensure even number (ch0 and ch1 always even length)
                start = (len(y) - ny) // 2
                y = y[int(start):int(start + ny)]

            if self.frame == "OD" or self.frame == "Column Density":
                if self.fitK is not None:
                    if self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac 2D"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac fixed betamu"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Azimuthal Average Gauss"
                    ):
                        self.figs.plotSliceUpdate(
                            x, [Sx, Fx], np.arange(len(R)), [R, RG, RF]
                        )
                    elif self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Gaussian Mask Sigma"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Gaussian Mask Sigma No Rot"
                    ) or self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Twisted Gauss (Mask) Int"
                    ):
                        self.figs.plotSliceUpdate(x, [Sx, Fx, Fxa], y, [Sy, Fy, Fya])
                    elif self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Fermi-Dirac 2D Int"
                    ):
                        # 2nd plot is x instead of y: it's for integration along y so you want to plot the x axis
                        self.figs.plotSliceUpdate(x, [Sx, Fx, Fy], x, [Sy])
                    
                    elif self.fitK.fitFunction == FIT_FUNCTIONS.index(
                        "Gauss (Mask) Int"
                    ):
                        # 2nd plot is x instead of y: it's for integration along y so you want to plot the x axis
                        self.figs.plotSliceUpdate(x, [Sx, Fx, Fxa], x, [Sy])
                    
                    else:
                        self.figs.plotSliceUpdate(x, [Sx, Fx], y, [Sy, Fy])

        if self.figs.plotTools.rbSelect.isChecked():

            box = self.pRb if self.av.b3_defr.isChecked() else None
            
            x = self.odRb.xRange0
            y = self.odRb.xRange1

            try:
                ch0 = self.fitRb.slices.ch0
                ch1 = self.fitRb.slices.ch1
                Sx = self.fitRb.slices.points0
                Sy = self.fitRb.slices.points1
                Fx = self.fitRb.slices.fit0
                Fxa = self.fitRb.slices.fit0a
                Fy = self.fitRb.slices.fit1
                Fya = self.fitRb.slices.fit1a
                
                # Initialized as 'None' in case the fit function does not have radial slices
                R = self.fitRb.slices.radSlice
                RG = self.fitRb.slices.radSliceFitGauss
                RF = self.fitRb.slices.radSliceFit

                if hasattr(self.fitRb.slices, "fit0Gauss"):
                    FxGauss = self.fitRb.slices.fit0Gauss
                else:
                    FxGauss = None

                if hasattr(self.fitRb.slices, "fit1Gauss"):
                    FyGauss = self.fitRb.slices.fit1Gauss
                else:
                    FyGauss = None

                if box is None: # old box usage below, not sure what the old usage of the box was
                    if hasattr(self.fitRb, "box"):
                        box = self.fitRb.box
                    else:
                        box = None
            except Exception as e:
                ch0 = None
                ch1 = None
                Sx = None
                Sy = None
                Fx = None
                Fxa = None
                Fy = None
                Fya = None
                FxGauss = None
                FyGauss = None
                # box = None
                self.figs.ax1.cla()
                self.figs.ax2.cla()
                print(e)

            if self.frame == "OD":
                image = self.odRb.ODCorrected
            elif self.frame == "Column Density":
                image = self.odRb.n
            else:
                image = frames[species[1]][self.frame][y[0] : y[-1], x[0] : x[-1]]
            self.figs.plotUpdate(x, y, image, ch0, ch1, box)

            if self.frame == "OD" or self.frame == "Column Density":
                if FxGauss is not None and FyGauss is not None:
                    self.figs.plotSliceUpdate(
                        x, [Sx, Fx, FxGauss], y, [Sy, Fy, FyGauss]
                    )
                else:
                    # This part will give an error if the fit isn't defined
                    # it is probably better to nest it in the next if statement
                    # but I opted to keep it like this since it does not break Imfit - Tim, 3/21/26
                    self.figs.plotSliceUpdate(x, [Sx, Fx], y, [Sy, Fy])
                if self.fitRb is not None:
                    if self.fitRb.fitFunction == FIT_FUNCTIONS.index(
                        "Gaussian Mask Sigma"
                    ) or self.fitRb.fitFunction == FIT_FUNCTIONS.index(
                        "Gaussian Mask Sigma No Rot"
                    ):
                        self.figs.plotSliceUpdate(x, [Sx, Fx, Fxa], y, [Sy, Fy, Fya])
                    elif self.fitRb.fitFunction == FIT_FUNCTIONS.index(
                        "Azimuthal Average Gauss"
                    ):
                        self.figs.plotSliceUpdate(
                            x, [Sx, Fx], np.arange(len(R)), [R, RG, RF]
                        )

    def passCamToROI(self):
        self.roi.setDefaultRegion(self.mode)
        self.roi.setCsat(self.mode)

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
        v0.addStretch(2)
        v0.addWidget(gb1)
        v0.addStretch(1)

        gb2 = QtWidgets.QGroupBox("Region Selection and Csat")
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

        gb4 = QtWidgets.QGroupBox("Averaging and Defringing")
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
        # # Added for Darkmode style:
        # p.setColor(QtGui.QPalette.Window, QtGui.QColor(53, 53, 53))
        # p.setColor(QtGui.QPalette.WindowText, QtCore.Qt.white)
        # p.setColor(QtGui.QPalette.Base, QtGui.QColor(25, 25, 25))
        # p.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(53, 53, 53))
        # p.setColor(QtGui.QPalette.ToolTipBase, QtCore.Qt.black)
        # p.setColor(QtGui.QPalette.ToolTipText, QtCore.Qt.white)
        # p.setColor(QtGui.QPalette.Text, QtCore.Qt.white)
        # p.setColor(QtGui.QPalette.Button, QtGui.QColor(53, 53, 53))
        # p.setColor(QtGui.QPalette.ButtonText, QtCore.Qt.white)
        # p.setColor(QtGui.QPalette.BrightText, QtCore.Qt.red)
        # p.setColor(QtGui.QPalette.Link, QtGui.QColor(42, 130, 218))
        # p.setColor(QtGui.QPalette.Highlight, QtGui.QColor(42, 130, 218))
        # p.setColor(QtGui.QPalette.HighlightedText, QtCore.Qt.black)
        # p.setColor(self.mainWidget.backgroundRole(), QtGui.QColor('#33373B'))  #QtCore.Qt.black)
        # p.setColor(self.mainWidget.backgroundRole(), QtCore.Qt.white)  #)
        p.setColor(self.mainWidget.backgroundRole(), QtGui.QColor("#fcfbfd"))
        self.mainWidget.setStyleSheet(self.getStyleSheet("./lib/styles.qss"))
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
                f = open(xp[0], "w")
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
    app.setStyle("Fusion")
    w = imfitDue()
    w.setGeometry(100, 100, 1200, 600)

    try:
        font = QtGui.QFont("Arial", 8)
        app.setFont(font)
    except Exception as e:
        raise (e)

    appico = QtGui.QIcon("IconKRb.ico")
    w.setWindowIcon(appico)

    w.show()
    sys.exit(app.exec_())

import numpy as np
import copy
from imfitHelpers import *
from imfitDefaults import *
from imfitFunctions import *

from scipy.optimize import least_squares 
from scipy.interpolate import interp2d

class calcOD():

    def __init__(self,data,atom,region=[]):

        ### Note that the region is passed as [x0, y0, xcrop, ycrop]
        ### Crop is symmetric about center (x0,y0)


        self.data = data
        self.xRange0 = None
        self.xRange1 = None

    
        if checkAtom(atom) < 0:
            pass
        else:
            self.atom = checkAtom(atom)
            if len(region) == 4:
                self.xCenter0 = region[0]
                self.xCenter1 = region[1]
                self.xCrop0 = region[2]
                self.xCrop1 = region[3]
                self.updateAll()
            else:
                print('The region of interest was not properly defined.')


    def setRegion(self, region):
        if len(region) == 4:
            self.xCenter0 = region[0]
            self.xCenter1 = region[1]
            self.xCrop0 = region[2]
            self.xCrop1 = region[3]
            self.updateAll()
        else:
            print('The region of interest was not properly defined.')

    
    
    def updateAll(self):
        self.defineROI()
        self.calculateOD()
        pass

    def calculateOD(self):
        with np.errstate(divide='ignore',invalid='ignore'):
            shadow = self.data.getFrame(self.atom,'Shadow')
            light = self.data.getFrame(self.atom,'Light')
            dark = self.data.getFrame(self.atom,'Dark')
    
            shadowCrop = cropArray(shadow, self.xRange1, self.xRange0)
            lightCrop = cropArray(light, self.xRange1, self.xRange0)
            darkCrop = cropArray(dark, self.xRange1, self.xRange0)
    
            self.OD = np.log((lightCrop-darkCrop)/(shadowCrop-darkCrop))
    
            if self.data.camera == CAMERA_NAME_PI:
                print(CAMERA_NAME_PI)
                isat = (self.data.bin+1.0)**2.0*ISAT_FLUX_PI[self.atom]*TPROBE_PI[self.atom]
                odsat = ODSAT_PI[self.atom]
            elif self.data.camera == CAMERA_NAME_IXON:
                print(CAMERA_NAME_IXON)
                isat = (self.data.bin+1.0)**2.0*ISAT_FLUX_IXON[self.atom]*TPROBE_IXON[self.atom]
                odsat = ODSAT_IXON[self.atom]
    
            ODmod = np.log((1.0-np.exp(-odsat))/(np.exp(-self.OD)-np.exp(-odsat)))
            self.ODCorrected = ODmod + (1-np.exp(-ODmod))*lightCrop/isat

            #Set all nans and infs to zero
            self.OD[np.isnan(self.OD)] = 0
            self.OD[np.isinf(self.OD)] = 0
            if self.data.camera == CAMERA_NAME_IXON:
                self.ODCorrected[np.isnan(self.ODCorrected)] = ODSAT_IXON[self.atom]
                self.ODCorrected[np.isinf(self.ODCorrected)] = ODSAT_IXON[self.atom]
            elif self.data.camera == CAMERA_NAME_PI:
                self.ODCorrected[np.isnan(self.ODCorrected)] = ODSAT_PI[self.atom]
                self.ODCorrected[np.isinf(self.ODCorrected)] = ODSAT_PI[self.atom]

    def defineROI(self):
        if self.xCenter0 > self.data.hImgSize or self.xCenter0 < 0:
            print('The horizontal center is out of range.')
            return -1
        elif self.xCenter1 > self.data.vImgSize or self.xCenter1 < 0:
            print('The vertical center is out of range.')
            return -1

        r0 = int(self.xCenter0 - np.floor(self.xCrop0/2))
        r1 = int(self.xCenter0 + np.floor(self.xCrop0/2))
        r2 = int(self.xCenter1 - np.floor(self.xCrop1/2))
        r3 = int(self.xCenter1 + np.floor(self.xCrop1/2))

        if r0 < 0:
            r0 = 0

        if r1 > self.data.hImgSize:
            r1 = self.data.hImgSize 

        if r2 < 0:
            r2 = 0

        if r3 > self.data.vImgSize:
            r3 = self.data.vImgSize

        self.xRange0 = range(r0,r1)
        self.xRange1 = range(r2,r3)
        
class fitOD():

    def __init__(self, odImage, fitFunction, atom, TOF, pxl):

        self.atom = atom
        self.TOF = TOF
        self.pxl = pxl*(1.0 + odImage.data.bin)

        self.fitFunction = None
        self.odImage = odImage

        self.fitData = None
        self.fitDataConf = None
        self.fittedImage = None


        self.slices = self.fitSlices()
        
        self.setFitFunction(fitFunction)
        self.fitODImage()

    class fitSlices():
        def __init__(self):
            self.points0 = []
            self.fit0 = []
            self.ch0 = None
            self.points1 = []
            self.fit1 = []
            self.ch1 = None
            self.radSlice = None
            self.radSliceFit = None
            self.radSliceFitGauss = None

    def setFitFunction(self, fitFunction):
       
        if not isinstance(fitFunction, basestring):
            print('Fit function must be one of: ' + ', '.join(FIT_FUNCTIONS)) 
            return -1

        if fitFunction not in FIT_FUNCTIONS:
            print('Fit function must be one of: ' + ', '.join(FIT_FUNCTIONS)) 
            return -1

        else:
            self.fitFunction = FIT_FUNCTIONS.index(fitFunction)

    def fitODImage(self):

        I0,I1 = np.unravel_index(self.odImage.ODCorrected.argmax(),self.odImage.ODCorrected.shape)
        M = self.odImage.ODCorrected[I0,I1]
        M *= (M>0) #Make sure M isn't out of bounds 

        if self.fitFunction == FIT_FUNCTIONS.index('Rotated Gaussian'):
            
            # Gaussian fit

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
    
            
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0]
            pUpper = [np.inf, 15.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1]), np.pi/2.0]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -np.pi/2.0]
            p0 = checkGuess(p0,pUpper,pLower)

            resLSQ = least_squares(gaussian, p0, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
            

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussian(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)

            
            ### Get radial average
            
            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)
            

            ####### Calculate slices through fit #####

            f = interp2d(self.odImage.xRange0, self.odImage.xRange1, self.odImage.ODCorrected, kind='cubic')


            m0 = np.tan(np.pi/2.0 - resLSQ.x[6])
            m1 = -np.tan(resLSQ.x[6])

            if abs(m0) > abs(m1):
                # Ensure the that lower slope is always along x
                m0, m1 = m1, m0


            b0 = -m0*resLSQ.x[2] + resLSQ.x[4]
            ch0 = np.asarray(self.odImage.xRange0)*m0 + b0
            b1 = -m1*resLSQ.x[2] + resLSQ.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1)/m1 


            for k in range(len(self.odImage.xRange0)):
                self.slices.points0.append(f(self.odImage.xRange0[k],ch0[k])[0])
                self.slices.fit0.append(gaussian(resLSQ.x,[r[0][k], ch0[k]],0.0)[0])
            self.slices.ch0 = ch0
 
            for k in range(len(self.odImage.xRange1)):
                self.slices.points1.append(f(ch1[k],self.odImage.xRange1[k])[0])
                self.slices.fit1.append(gaussian(resLSQ.x,[ch1[k], r[1][k]],0)[0])
            self.slices.ch1 = ch1

        elif self.fitFunction == FIT_FUNCTIONS.index('Gaussian'):

            # Gaussian fit without rotation

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
    
            
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]
            pUpper = [np.inf, 15.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0,pUpper,pLower)

            resLSQ = least_squares(gaussianNoRot, p0, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
            

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussianNoRot(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)

            
            ### Get radial average
            
            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1,:]
            self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1,:]

            self.slices.points1 = self.odImage.ODCorrected[:,I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:,I0]


        elif self.fitFunction == FIT_FUNCTIONS.index('Bigaussian'):
            #Bi Gaussian fit
            ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
    
            ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]
            p0 = [0, M/2.0, 20.0, 20.0, M/2.0, 8.0, 8.0, self.odImage.xRange0[I1], self.odImage.xRange1[I0]]
            pUpper = [np.inf, 8.0, len(r[0]), len(r[1]), 8.0, len(r[0]), len(r[1]), np.max(r[0]), np.max(r[1])]
            pLower = [-np.inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            p0 = checkGuess(p0,pUpper,pLower)

            resLSQ = least_squares(doubleGaussian, p0, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
            #self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = doubleGaussian(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)
                
            
            ### Make sure the BEC comes first in the list

            if self.fitData[2] + self.fitData[3] > self.fitData[5] + self.fitData[6]:
                t = copy.deepcopy(self.fitData[4:7])

                self.fitData[4:7] = self.fitData[1:4]
                self.fitData[1:4] = t
            else:
                pass

            ### Calculate radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[7]))
            I1 = self.odImage.xRange1.index(int(self.fitData[8]))

            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected,[I0, I1])
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, [I0, I1])
            self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1,:]
            self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1,:]

            self.slices.points1 = self.odImage.ODCorrected[:,I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:,I0]


        elif self.fitFunction == FIT_FUNCTIONS.index('Fermi-Dirac'):
            
            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            
            
            #Fermi--Dirac fit
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]            
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0] #An initial q of 0 corresponds to T/TF=0.56
            pUpper = [np.inf, 900.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1]), np.inf]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -np.inf]
            p0 = checkGuess(p0,pUpper,pLower)

            resLSQ = least_squares(fermiDirac, p0, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = fermiDirac(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)

            #Gaussian fit
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0] 

            print(p0)
            pUpper = [np.inf, 15.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1]), 0.00001]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -0.00001]

            resLSQ = least_squares(gaussian, p0, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussian(self.fitDataGauss, r, 0).reshape(self.odImage.ODCorrected.shape)

            #######################################################################################################################

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            # azAverage -- I'm not sure what the correct index should be for azAverage... 
            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = azimuthalAverage(self.fittedImageGauss, center)

            ### Calculate slices through fit (No rotation)
            self.slices.points0 = self.odImage.ODCorrected[I1,:]
            self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1,:] 

            self.slices.points1 = self.odImage.ODCorrected[:,I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
            self.slices.fit1 =  self.fittedImage[:,I0]

        elif self.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):

            # Band mapping function in the vertical direction

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1


            N0 = len(self.odImage.xRange0)
            N1 = len(self.odImage.xRange1)

            ### Parameters: [Offset, A0, A1, A2, wy, yc, wx, xc]
            p0 = [0, M, M*0.05, 0, 10.0, self.odImage.xRange1[N1/2], 40, self.odImage.xRange0[N0/2]]
            pUpper = [np.inf, 15.0, 15.0, 15.0, len(r[1]), np.max(r[1]), len(r[0]), np.max(r[0])]
            pLower = [-np.inf, 0.0, 0.0, 0.0,  0.0, np.min(r[1]), 0.0, np.min(r[0])]
            p0 = checkGuess(p0,pUpper,pLower)

            imageDetails = [self.atom, self.TOF, self.pxl]

            resLSQ = least_squares(bandmapV, p0, args=(r,self.odImage.ODCorrected, imageDetails),bounds=(pLower,pUpper))
            

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = bandmapV(resLSQ.x, r, 0,imageDetails).reshape(self.odImage.ODCorrected.shape)            

            ### Get radial average
            
            I0 = self.odImage.xRange0.index(int(self.fitData[7]))
            I1 = self.odImage.xRange1.index(int(self.fitData[5]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1,:]
            self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1,:]

            self.slices.points1 = self.odImage.ODCorrected[:,I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:,I0]
        
        else:
            print('Fit function undefined! Something went wrong!')
            return -1

class processFitResult():

    def __init__(self,fitObject, imagePath='Axial'):

        self.fitObject = fitObject 
        self.bin = self.fitObject.odImage.data.bin
        self.atom = self.fitObject.odImage.atom
        self.imagePath = IMAGING_PATHS.index(imagePath)
        self.pixelSize = PIXEL_SIZES[self.imagePath]

        self.data = None

        self.getResults()

    def getResults(self):

        
        if self.fitObject.fitFunction == FIT_FUNCTIONS.index('Rotated Gaussian'):
            
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*(self.bin+1.0)*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*(self.bin+1.0)*self.pixelSize,
                    'angle' : self.fitObject.fitData[6]*180.0/3.141592
                    }

            self.data = ['fileName', r['peakOD'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']]

            rErr = {
                    'offset' : self.fitObject.fitDataConf[0],
                    'peakOD' : self.fitObject.fitDataConf[1],
                    'x0' : self.fitObject.fitDataConf[2],
                    'y0' : self.fitObject.fitDataConf[4],
                    'wx' : self.fitObject.fitDataConf[3]*(self.bin+1.0)*self.pixelSize,
                    'wy' : self.fitObject.fitDataConf[5]*(self.bin+1.0)*self.pixelSize,
                    'angle' : self.fitObject.fitDataConf[6]*180.0/3.141592
                    }

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Gaussian'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*(self.bin+1.0)*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*(self.bin+1.0)*self.pixelSize,
                    'angle' : 0
                    }

            self.data = ['fileName', r['peakOD'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']]            

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Bigaussian'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakODBEC' : self.fitObject.fitData[1],
                    'wxBEC' : self.fitObject.fitData[2]*(self.bin+1.0)*self.pixelSize,
                    'wyBEC' : self.fitObject.fitData[3]*(self.bin+1.0)*self.pixelSize,
                    'peakODThermal' : self.fitObject.fitData[4],
                    'wxThermal' : self.fitObject.fitData[5]*(self.bin+1.0)*self.pixelSize,
                    'wyThermal' : self.fitObject.fitData[6]*(self.bin+1.0)*self.pixelSize,
                    'x0': self.fitObject.fitData[7],
                    'y0': self.fitObject.fitData[8],
                    }

            self.data = ['fileName', r['peakODBEC'], r['wxBEC'], r['wyBEC'], r['peakODThermal'], r['wxThermal'], r['wyThermal'], r['x0'], r['y0'], r['offset']]

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Fermi-Dirac'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*(self.bin+1.0)*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*(self.bin+1.0)*self.pixelSize,
                    'q' : self.fitObject.fitData[6],
                    'TTF': getTTF(self.fitObject)[0][0],
                    'wxClassical' : self.fitObject.fitDataGauss[3]*(self.bin+1.0)*self.pixelSize,
                    'wyClassical' : self.fitObject.fitDataGauss[5]*(self.bin+1.0)*self.pixelSize,
                    'peakODClassical' : self.fitObject.fitDataGauss[1]
                    }

            self.data = ['fileName', r['peakODClassical'], r['wxClassical'], r['wyClassical'], r['peakOD'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['TTF']]

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'Band0' : self.fitObject.fitData[1],
                    'Band1' : self.fitObject.fitData[2],
                    'Band2' : self.fitObject.fitData[3],
                    'wy' : self.fitObject.fitData[4]*(self.bin+1.0)*self.pixelSize,
                    'wx' : self.fitObject.fitData[6]*(self.bin+1.0)*self.pixelSize,
                    'x0': self.fitObject.fitData[7],
                    'y0': self.fitObject.fitData[5],
                    'TOF': self.fitObject.TOF,
                    }

            self.data = ['fileName', 'species', r['Band0'], r['Band1'], r['Band2'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['TOF']]

        else:
            print('Fit function undefined! Something went wrong!')
            return -1

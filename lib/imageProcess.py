from math import isnan

from scipy.ndimage.measurements import histogram, maximum
from lib.imfitDefaults import DEFAULT_MODE, IMFIT_MODES
import numpy as np
import copy
from lib.imfitHelpers import *
from lib.imfitDefaults import *
from lib.imfitFunctions import *

from scipy.optimize import least_squares 
from scipy.interpolate import interp2d

from skimage.feature import peak_local_max
from scipy import ndimage as ndimage

class calcOD():

    def __init__(self,data,species,mode,region=[]):

        ### Note that the region is passed as [x0, y0, xcrop, ycrop]
        ### Crop is symmetric about center (x0,y0)

        self.data = data
        self.xRange0 = None
        self.xRange1 = None

        self.mode = mode
        self.config = IMFIT_MODES[mode]
        self.species = species
    
        if len(region) == 4:
            self.xCenter0 = region[0]
            self.xCenter1 = region[1]
            self.xCrop0 = region[2]
            self.xCrop1 = region[3]
            self.updateAll()
        else:
            print('The region of interest was not properly defined.')


        # if checkAtom(species) < 0:
        #     pass
        # else:
        #     self.atom = checkAtom(species)
        #     if len(region) == 4:
        #         self.xCenter0 = region[0]
        #         self.xCenter1 = region[1]
        #         self.xCrop0 = region[2]
        #         self.xCrop1 = region[3]
        #         self.updateAll()
        #     else:
        #         print('The region of interest was not properly defined.')


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
            shadow = self.data.getFrame(self.species,'Shadow')
            light = self.data.getFrame(self.species,'Light')
            dark = self.data.getFrame(self.species,'Dark')
    
            shadowCrop = cropArray(shadow, self.xRange1, self.xRange0)
            lightCrop = cropArray(light, self.xRange1, self.xRange0)
            darkCrop = cropArray(dark, self.xRange1, self.xRange0)

            s1 = shadowCrop - darkCrop
            s2 = lightCrop - darkCrop
            self.OD = -np.log(s1/s2)

            Ceff = self.config['CSat'][self.species]
            bins = self.data.bin
            Ceff *= float(bins**2)

            self.ODCorrected = self.OD + (s2 - s1)/Ceff

            #Set all nans and infs to zero
            self.ODCorrected[np.isnan(self.ODCorrected)] = 0
            self.ODCorrected[np.isinf(self.ODCorrected)] = 0

            # Calculate the column density, assuming zero detuning; see Pappa et al, NJP (2011)
            delta = 0 # Detuning; assumed to be zero
            Omega = 4*np.pi*np.sin(np.arcsin(self.config["NA"])/2)**2 # Fraction of flourescence collected
            sigma0 = SIGMA_0[self.species] # Resonant cross section at I/Isat = 0, in um^2
            s0 = s2/(bins**2 * self.config['CSat'][self.species])
            Tabs = s1/s2
            self.n = (1 + delta**2)/((1 - Omega) * sigma0) * (- np.log(Tabs) + s0/(1 + delta**2) * (1 - Tabs))
            self.n /= EFF[self.species]
            self.n[np.isnan(self.n)] = 0
            self.n[np.isinf(self.n)] = 0

            self.nerr = np.sqrt((s2+s1)/(s2**3 * s1)) * (s2*(1 + delta**2) + s0*s1) / (sigma0 * (1 - Omega))
            self.nerr[np.isnan(self.nerr)] = 0
            self.nerr[np.isinf(self.nerr)] = 0
  
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

    def __init__(self, mode, odImage, fitFunction, species, TOF, pxl):

        self.species = species
        self.TOF = TOF

        self.pxl = pxl*(1.0 + odImage.data.bin)

        self.fitFunction = None
        self.odImage = odImage

        self.fitData = None
        self.fitDataConf = None
        self.fittedImage = None

        self.slices = self.fitSlices()
        
        self.mode = mode
        self.config = IMFIT_MODES[mode]

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
       
        if not isinstance(fitFunction, str):
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

        if M > 10:
            M = 10

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


        elif self.fitFunction == FIT_FUNCTIONS.index('Twisted Gaussian'):
            
            # Gaussian fit, rotated by VERT_TRAP_ANGLE

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
    
            
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]
            pUpper = [np.inf, 15.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0,pUpper,pLower)

            angle = IMFIT_MODES[self.mode]['Fit angle']
            resLSQ = least_squares(gaussianNoRotTwist, p0, args=(r,self.odImage.ODCorrected,angle),bounds=(pLower,pUpper))

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussianNoRotTwist(resLSQ.x, r, 0, angle).reshape(self.odImage.ODCorrected.shape)

            
            ### Get radial average
            
            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)
            

            ####### Calculate slices through fit #####

            f = interp2d(self.odImage.xRange0, self.odImage.xRange1, self.odImage.ODCorrected, kind='cubic')


            m0 = np.tan(np.pi/2.0 - angle*np.pi/180.)
            m1 = -np.tan(angle*np.pi/180.)

            if abs(m0) > abs(m1):
                # Ensure the that lower slope is always along x
                m0, m1 = m1, m0


            b0 = -m0*resLSQ.x[2] + resLSQ.x[4]
            ch0 = np.asarray(self.odImage.xRange0)*m0 + b0
            b1 = -m1*resLSQ.x[2] + resLSQ.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1)/m1 


            for k in range(len(self.odImage.xRange0)):
                self.slices.points0.append(f(self.odImage.xRange0[k],ch0[k])[0])
                self.slices.fit0.append(gaussianNoRotTwist(resLSQ.x,[r[0][k], ch0[k]],0.0,angle)[0])
            self.slices.ch0 = ch0
 
            for k in range(len(self.odImage.xRange1)):
                self.slices.points1.append(f(ch1[k],self.odImage.xRange1[k])[0])
                self.slices.fit1.append(gaussianNoRotTwist(resLSQ.x,[ch1[k], r[1][k]],0,angle)[0])
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

        elif self.fitFunction == FIT_FUNCTIONS.index('Gaussian w/ Gradient'):

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D,Y2D = np.meshgrid(np.arange(dx),np.arange(dy))
                A = np.matrix(np.column_stack((X2D.ravel(),Y2D.ravel(),np.ones(dx*dy))))
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg=np.reshape(C*A.T, (dy,dx))
                return np.asarray(od - bg)

            # Gaussian fit without rotation
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg,5,mode='constant')
    
            
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0, 0]
            pUpper = [np.inf, 15.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1]), 40, 40]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -40, -40]
            p0 = checkGuess(p0,pUpper,pLower)

            pks=peak_local_max(blur, min_distance=20,exclude_border=2, num_peaks=3)
            guesses = [p0]
            for pk in pks:
                yc = pk[0]
                xc = pk[1]
                peak = data[yc, xc]
                offset = np.mean(data)
                if peak > 0:
                    (sigx, sigy) = 15, 15
                    guess = [offset, peak, xmin+xc, sigx, ymin+yc, sigy, 0.0, 0.0]
                    guesses.append(checkGuess(guess,pUpper,pLower))

            best_fit = None
            best_guess = np.inf		
            for guess in guesses:
                try:
                    # print("Trying guess: {}".format(guess))
                    res = least_squares(gaussianNoRotGradient, guess, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
                    # print("Cost: {}".format(res.cost))
                    if not res.success:
                        print("Warning: fit did not converge.")
                    elif res.cost < best_guess:
                        best_guess = res.cost
                        best_fit = res
                except ValueError as e:
                    print("Error fitting image: {}".format(e))

            resLSQ = best_fit
            # print("Best fit: {}".format(resLSQ))
            if resLSQ is not None:
                self.fitDataConf = confidenceIntervals(resLSQ)
                self.fitData = resLSQ.x
                self.fittedImage = gaussianNoRotGradient(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)

            
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

            print("Done with fit function!")
        
        elif self.fitFunction == FIT_FUNCTIONS.index('Gaussian Fixed'):

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D,Y2D = np.meshgrid(np.arange(dx),np.arange(dy))
                A = np.matrix(np.column_stack((X2D.ravel(),Y2D.ravel(),np.ones(dx*dy))))
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg=np.reshape(C*A.T, (dy,dx))
                return np.asarray(od - bg)

            # Gaussian fit without rotation
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]

            r = [None,None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])
            xmax = np.max(r[0])
            ymax = np.max(r[1])
            xc = 0.5*(xmin + xmax)
            yc = 0.5*(ymin + ymax)

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg,5,mode='constant')
    
            
            p0 = [0, M, xc, 15, yc, 3, 0, 0]
            pUpper = [np.inf, 15.0, xc + 2, 40, yc + 2, 6, 40, 40]
            pLower = [-np.inf, 0.0, xc - 2, 4, yc - 2, 0.01, -40, -40]
            p0 = checkGuess(p0,pUpper,pLower)

            guesses = [p0]

            best_fit = None
            best_guess = np.inf		
            for guess in guesses:
                try:
                    print("Trying guess: {}".format(guess))
                    res = least_squares(gaussianNoRotGradient, guess, args=(r,self.odImage.ODCorrected),bounds=(pLower,pUpper))
                    # print("Cost: {}".format(res.cost))
                    if not res.success:
                        print("Warning: fit did not converge.")
                    elif res.cost < best_guess:
                        best_guess = res.cost
                        best_fit = res
                except ValueError as e:
                    print("Error fitting image: {}".format(e))

            resLSQ = best_fit
            # print("Best fit: {}".format(resLSQ))
            if resLSQ is not None:
                self.fitDataConf = confidenceIntervals(resLSQ)
                self.fitData = resLSQ.x
                self.fittedImage = gaussianNoRotGradient(resLSQ.x, r, 0).reshape(self.odImage.ODCorrected.shape)

            
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

            print("Done with fit function!")


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
            x_coerced = min(max(self.fitData[7], self.odImage.xRange0.start), self.odImage.xRange0.stop-1)
            y_coerced = min(max(self.fitData[8], self.odImage.xRange1.start), self.odImage.xRange1.stop-1)
            
            I0 = self.odImage.xRange0.index(int(x_coerced))
            I1 = self.odImage.xRange1.index(int(y_coerced))

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

        # elif self.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):
        #     # TODO: Make this depend properly on the species. Disabled for now.
        #     # Band mapping function in the vertical direction

        #     r = [None,None]
        #     r[0] = self.odImage.xRange0
        #     r[1] = self.odImage.xRange1


        #     N0 = len(self.odImage.xRange0)
        #     N1 = len(self.odImage.xRange1)

        #     ### Parameters: [Offset, A0, A1, A2, wy, yc, wx, xc]
        #     p0 = [0, M, M*0.05, 0, 10.0, self.odImage.xRange1[N1/2], 40, self.odImage.xRange0[N0/2]]
        #     pUpper = [np.inf, 15.0, 15.0, 15.0, len(r[1]), np.max(r[1]), len(r[0]), np.max(r[0])]
        #     pLower = [-np.inf, 0.0, 0.0, 0.0,  0.0, np.min(r[1]), 0.0, np.min(r[0])]
        #     p0 = checkGuess(p0,pUpper,pLower)

        #     imageDetails = [self.species, self.TOF, self.pxl]

        #     resLSQ = least_squares(bandmapV, p0, args=(r,self.odImage.ODCorrected, imageDetails),bounds=(pLower,pUpper))
            

        #     self.fitDataConf = confidenceIntervals(resLSQ)
        #     self.fitData = resLSQ.x
        #     self.fittedImage = bandmapV(resLSQ.x, r, 0,imageDetails).reshape(self.odImage.ODCorrected.shape)            

        #     ### Get radial average
            
        #     I0 = self.odImage.xRange0.index(int(self.fitData[7]))
        #     I1 = self.odImage.xRange1.index(int(self.fitData[5]))

        #     center = [I0, I1]
        #     self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
        #     self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
        #     self.slices.radSliceFitGauss = [0]*len(self.slices.radSlice)

        #     ### Calculate slices through fit

        #     self.slices.points0 = self.odImage.ODCorrected[I1,:]
        #     self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
        #     self.slices.fit0 = self.fittedImage[I1,:]

        #     self.slices.points1 = self.odImage.ODCorrected[:,I0]
        #     self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
        #     self.slices.fit1 = self.fittedImage[:,I0]

        elif self.fitFunction == FIT_FUNCTIONS.index('Integrate'):
            # Integration of number from computed column density

            # Calculate average number density in border and subtract from rest of image
            border = max(1, int(min(min(self.odImage.n.shape)/10, 5)))
            print("Border: {}".format(border))
            interior = self.odImage.n[border:-border, border:-border]
            offset = (self.odImage.n.sum() - interior.sum()) / (self.odImage.n.size - interior.size)
            print("Offset: {}".format(offset))
            interior -= offset

            interior_err = self.odImage.nerr[border:-border, border:-border]

            # Compute the number by summing the pixels and multiplying by the pixel area

            raw_number = interior.sum()
            number = raw_number * (self.config['Pixel Size'] * self.odImage.data.bin)**2
            print("Number: {}".format(number))

            raw_error = np.sqrt((interior_err**2).sum())
            error = raw_error * (self.config['Pixel Size'] * self.odImage.data.bin)**2
            print("Error: {}".format(error))

            # Compute the central position and size
            xrange = np.arange(interior.shape[1])
            yrange = np.arange(interior.shape[0])
            xc = np.sum(interior.sum(axis=0) * xrange) / raw_number
            yc = np.sum(interior.sum(axis=1) * yrange) / raw_number

            # These calculations of moments aren't exact, since the contents of the square root can be negative since the column density can be negative in some pixels.
            # Therefore, all pixels less than 10 percent of the peak are dropped.
            sumx = interior.sum(axis=0)
            sumy = interior.sum(axis=1)
            minx = maximum(sumx)/10
            miny = maximum(sumy)/10
            sigx = np.sqrt(np.sum((sumx * (xrange - xc)**2)[sumx > minx]) / np.sum(sumx[sumx > minx]))
            sigy = np.sqrt(np.sum((sumy * (yrange - yc)**2)[sumy > miny]) / np.sum(sumy[sumy > miny]))

            xc += border + self.odImage.xRange0.start
            yc += border + self.odImage.xRange1.start

            xc = min(max(self.odImage.xRange0.start, xc),self.odImage.xRange0.stop-1)
            yc = min(max(self.odImage.xRange1.start, yc),self.odImage.xRange1.stop-1)

            print("xc: {} px".format(xc))
            print("yc: {} px".format(yc))
            print("sigx: {} px".format(sigx))
            print("sigy: {} px".format(sigy))

            if isnan(sigx):
                print("Sigma x invalid: setting to zero.")
                sigx = 0
            
            if isnan(sigy):
                print("Sigma y invalid: setting to zero.")
                sigy = 0

            self.fitData = [offset, number, error, xc, yc, sigx, sigy]
            self.fitDataConf = None
            self.fittedImage = None

            
            ### Get radial average
            
            I0 = self.odImage.xRange0.index(round(xc))
            I1 = self.odImage.xRange1.index(round(yc))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.n, center)
            self.slices.radSliceFit = None
            self.slices.radSliceFitGauss = None

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.n[I1,:]
            self.slices.ch0 = [self.odImage.xRange1[I1]]*len(self.odImage.xRange0)
            self.slices.fit0 = None

            self.slices.points1 = self.odImage.n[:,I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]]*len(self.odImage.xRange1)
            self.slices.fit1 = None
        
        else:
            print('Fit function undefined! Something went wrong!')
            return -1

class processFitResult():

    def __init__(self,fitObject, mode):

        self.fitObject = fitObject 
        self.bin = self.fitObject.odImage.data.bin
        # self.atom = self.fitObject.odImage.atom
        self.config = IMFIT_MODES[mode]
        self.pixelSize = self.config['Pixel Size']
        self.angle = self.config['Fit angle']

        self.data = None
        self.data_dict = None

        self.getResults()

    def getResults(self):        
        if self.fitObject.fitFunction == FIT_FUNCTIONS.index('Rotated Gaussian'):
            
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'angle' : self.fitObject.fitData[6]*180.0/np.pi,
                    'dODdx' : 0,
                    'dODdy' : 0,
                    }

            self.data =  ['fileName', r['peakOD'], r['dODdx'], r['dODdy'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']]
            self.data_dict = r

            rErr = {
                    'offset' : self.fitObject.fitDataConf[0],
                    'peakOD' : self.fitObject.fitDataConf[1],
                    'x0' : self.fitObject.fitDataConf[2],
                    'y0' : self.fitObject.fitDataConf[4],
                    'wx' : self.fitObject.fitDataConf[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitDataConf[5]*self.bin*self.pixelSize,
                    'angle' : self.fitObject.fitDataConf[6]*180.0/np.pi
                    }

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Gaussian'):
            
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'angle' : 0,
                    'dODdx' : 0,
                    'dODdy' : 0,
                    }

            self.data =['fileName', r['peakOD'], r['dODdx'], r['dODdy'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']] 
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Gaussian w/ Gradient') or self.fitObject.fitFunction == FIT_FUNCTIONS.index('Gaussian Fixed'):
            
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'dODdx' : self.fitObject.fitData[6]/self.bin*self.pixelSize,
                    'dODdy' : self.fitObject.fitData[7]/self.bin*self.pixelSize,
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'angle' : 0
                    }

            print(r)
            self.data = ['fileName', r['peakOD'], r['dODdx'], r['dODdy'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']]
            self.data_dict = r


        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Twisted Gaussian'):
            
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'angle' : self.angle,
                    'dODdx' : 0,
                    'dODdy' : 0,
                    }

            self.data = ['fileName', r['peakOD'], r['dODdx'], r['dODdy'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['angle']]        
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Bigaussian'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakODBEC' : self.fitObject.fitData[1],
                    'wxBEC' : self.fitObject.fitData[2]*self.bin*self.pixelSize,
                    'wyBEC' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'peakODThermal' : self.fitObject.fitData[4],
                    'wxThermal' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'wyThermal' : self.fitObject.fitData[6]*self.bin*self.pixelSize,
                    'x0': self.fitObject.fitData[7],
                    'y0': self.fitObject.fitData[8],
                    }

            self.data = ['fileName', r['peakODBEC'], r['wxBEC'], r['wyBEC'], r['peakODThermal'], r['wxThermal'], r['wyThermal'], r['x0'], r['y0'], r['offset']]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Fermi-Dirac'):

            r = {
                    'offset' : self.fitObject.fitData[0],
                    'peakOD' : self.fitObject.fitData[1],
                    'x0' : self.fitObject.fitData[2],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[3]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'q' : self.fitObject.fitData[6],
                    'TTF': getTTF(self.fitObject)[0][0],
                    'wxClassical' : self.fitObject.fitDataGauss[3]*self.bin*self.pixelSize,
                    'wyClassical' : self.fitObject.fitDataGauss[5]*self.bin*self.pixelSize,
                    'peakODClassical' : self.fitObject.fitDataGauss[1]
                    }

            self.data = ['fileName', r['peakODClassical'], r['wxClassical'], r['wyClassical'], r['peakOD'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['TTF']]
            self.data_dict = r

        # elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Vertical BandMap'):

        #     r = {
        #             'offset' : self.fitObject.fitData[0],
        #             'Band0' : self.fitObject.fitData[1],
        #             'Band1' : self.fitObject.fitData[2],
        #             'Band2' : self.fitObject.fitData[3],
        #             'wy' : self.fitObject.fitData[4]*self.bin*self.pixelSize,
        #             'wx' : self.fitObject.fitData[6]*self.bin*self.pixelSize,
        #             'x0': self.fitObject.fitData[7],
        #             'y0': self.fitObject.fitData[5],
        #             'TOF': self.fitObject.TOF,
        #             }

        #     self.data = ['fileName', 'species', r['Band0'], r['Band1'], r['Band2'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset'], r['TOF']]

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index('Integrate'):
            r = {
                    'offset' : self.fitObject.fitData[0],
                    'number' : self.fitObject.fitData[1],
                    'error' : self.fitObject.fitData[2],
                    'x0' : self.fitObject.fitData[3],
                    'y0' : self.fitObject.fitData[4],
                    'wx' : self.fitObject.fitData[5]*self.bin*self.pixelSize,
                    'wy' : self.fitObject.fitData[6]*self.bin*self.pixelSize,
                }
                
            self.data = ['fileName', r['number'], r['error'], r['wx'], r['wy'], r['x0'], r['y0'], r['offset']]
            self.data_dict = r

        else:
            print('Fit function undefined! Something went wrong!')
            return -1

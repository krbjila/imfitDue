from math import isnan

from scipy.ndimage.measurements import histogram, maximum
from lib.imfitDefaults import DEFAULT_MODE, IMFIT_MODES
import numpy as np
import copy
from lib.imfitHelpers import *
from lib.imfitDefaults import *
from lib.imfitFunctions import *
from lib.polylog import fermi_poly2

from scipy.optimize import least_squares
from scipy.interpolate import RectBivariateSpline

from skimage.feature import peak_local_max
from scipy import ndimage as ndimage


class calcOD:

    def __init__(self, data, species, mode, region=[], sigblur=0):

        ### Note that the region is passed as [x0, y0, xcrop, ycrop]
        ### Crop is symmetric about center (x0,y0)

        self.data = data
        self.sigblur = sigblur
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
            print("The region of interest was not properly defined.")

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
            print("The region of interest was not properly defined.")

    def updateAll(self):
        self.defineROI()
        self.calculateOD()
        pass

    def calculateOD(self):
        with np.errstate(divide="ignore", invalid="ignore"):
            shadow = self.data.getFrame(self.species, "Shadow")
            light = self.data.getFrame(self.species, "Light")
            dark = self.data.getFrame(self.species, "Dark")

            shadowCrop = cropArray(shadow, self.xRange1, self.xRange0)
            lightCrop = cropArray(light, self.xRange1, self.xRange0)
            darkCrop = cropArray(dark, self.xRange1, self.xRange0)

            s1 = shadowCrop - darkCrop
            s2 = lightCrop - darkCrop
            self.OD = -np.log(s1 / s2)

            # Correct OD for fluorescence and saturation
            # Detuning; assumed to be zero
            delta = 0
            # Fraction of fluorescence collected
            Omega = (1 - np.cos(np.arcsin(self.config["NA"]))) / 2
            # Resonant cross section at I/Isat = 0, in um^2
            sigma0 = SIGMA_0[self.species] * (
                2 if IMFIT_MODES[self.mode]["Image Path"] == "Vertical" else 1
            )  # times 2 for circular polarization

            Ceff = self.config["CSat"][self.species]
            bins = self.data.bin
            Ceff *= float(bins**2)

            self.ODCorrected = (self.OD + (s2 - s1) / Ceff) / (1 - Omega)

            # Set all nans and infs to zero
            self.ODCorrected[np.isnan(self.ODCorrected)] = 0
            self.ODCorrected[np.isinf(self.ODCorrected)] = 0

            if self.sigblur > 0:
                self.GaussBlur(self.sigblur)
            
            # uncomment next line if you want to save the corrected OD for debugging purposes
            # np.savetxt('output_array_2026-02-27.txt', self.ODCorrected, delimiter=',')

            # Calculate the column density, assuming zero detuning; see Pappa et al, NJP (2011)
            s0 = s2 / (bins**2 * self.config["CSat"][self.species])
            Tabs = s1 / s2
            self.n = (-np.log(Tabs) + s0 * (1 - Tabs)) / ((1 - Omega) * sigma0)
            self.n /= EFF[self.species]
            self.n[np.isnan(self.n)] = 0
            self.n[np.isinf(self.n)] = 0

            self.nerr = (
                (s0 * Tabs + (1 + delta ^ 2))
                / (sigma0 * (1 - Omega))
                * np.sqrt((s1 + s2) / (s1 * s2))
            )
            self.nerr[np.isnan(self.nerr)] = 0
            self.nerr[np.isinf(self.nerr)] = 0

    def defineROI(self):
        if self.xCenter0 > self.data.hImgSize or self.xCenter0 < 0:
            print("The horizontal center is out of range.")
            return -1
        elif self.xCenter1 > self.data.vImgSize or self.xCenter1 < 0:
            print("The vertical center is out of range.")
            return -1

        r0 = int(self.xCenter0 - np.floor(self.xCrop0 / 2))
        r1 = int(self.xCenter0 + np.floor(self.xCrop0 / 2))
        r2 = int(self.xCenter1 - np.floor(self.xCrop1 / 2))
        r3 = int(self.xCenter1 + np.floor(self.xCrop1 / 2))

        if r0 < 0:
            r0 = 0

        if r1 > self.data.hImgSize:
            r1 = self.data.hImgSize

        if r2 < 0:
            r2 = 0

        if r3 > self.data.vImgSize:
            r3 = self.data.vImgSize

        self.xRange0 = range(r0, r1)
        self.xRange1 = range(r2, r3)
    
    def GaussBlur(self, sigblur):
        """
        Apply a low-pass FFT filter to a 2D image.
        
        Parameters:
            img (ndarray): 2D numpy array
            sigblur (float): Standard deviation of the Gaussian kernel
        """
        img = self.ODCorrected
        self.ODCorrected = ndimage.gaussian_filter(img, sigma=sigblur)



class fitOD:

    def __init__(self, mode, odImage, fitFunction, species, TOF, pxl, **kwargs):

        self.species = species
        self.TOF = TOF

        self.pxl = pxl * (1.0 + odImage.data.bin)

        self.fitFunction = None
        self.odImage = odImage

        self.fitData = None
        self.fitDataConf = None
        self.fittedImage = None

        self.slices = self.fitSlices()

        self.mode = mode
        self.config = IMFIT_MODES[mode]

        for k, val in kwargs.items():
            if k == "WingRad":
                self.WingRad = val
            if k == "fx":
                self.fx = val
            if k == "fy":
                self.fy = val
            if k == "fz":
                self.fz = val
            if k == "mbemu":
                self.mbemu = val
            if k == "Nscaler":
                self.Nscaler = val

        self.setFitFunction(fitFunction)
        self.fitODImage()

    class fitSlices:
        def __init__(self):
            self.points0 = []
            self.fit0 = []
            self.fit0a = []  # Added for Gaussian Wing Fitting
            self.ch0 = None
            self.points1 = []
            self.fit1 = []
            self.fit1a = []  # Added for Gaussian Wing Fitting
            self.ch1 = None
            self.radSlice = None
            self.radSliceFit = None
            self.radSliceFitGauss = None

    def setFitFunction(self, fitFunction):

        if not isinstance(fitFunction, str):
            print("Fit function must be one of: " + ", ".join(FIT_FUNCTIONS))
            return -1

        if fitFunction not in FIT_FUNCTIONS:
            print("Fit function must be one of: " + ", ".join(FIT_FUNCTIONS))
            return -1

        else:
            self.fitFunction = FIT_FUNCTIONS.index(fitFunction)

    def fitODImage(self):

        I0, I1 = np.unravel_index(
            self.odImage.ODCorrected.argmax(), self.odImage.ODCorrected.shape
        )
        M = self.odImage.ODCorrected[I0, I1]
        M *= M > 0  # Make sure M isn't out of bounds

        if M > 10:
            M = 10

        if self.fitFunction == FIT_FUNCTIONS.index("Rotated Gaussian"):

            # Gaussian fit

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0]
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                np.pi / 2.0,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -np.pi / 2.0]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                gaussian,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussian(resLSQ.x, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ####### Calculate slices through fit #####

            f = RectBivariateSpline(
                self.odImage.xRange1,
                self.odImage.xRange0,
                self.odImage.ODCorrected,
            )

            m0 = np.tan(np.pi / 2.0 - resLSQ.x[6])
            m1 = -np.tan(resLSQ.x[6])

            if abs(m0) > abs(m1):
                # Ensure the that lower slope is always along x
                m0, m1 = m1, m0

            b0 = -m0 * resLSQ.x[2] + resLSQ.x[4]
            ch0 = np.asarray(self.odImage.xRange0) * m0 + b0
            b1 = -m1 * resLSQ.x[2] + resLSQ.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1) / m1

            for k in range(len(self.odImage.xRange0)):
                self.slices.points0.append(f(ch0[k], self.odImage.xRange0[k])[0])
                self.slices.fit0.append(gaussian(resLSQ.x, [r[0][k], ch0[k]], 0.0)[0])
            self.slices.ch0 = ch0

            for k in range(len(self.odImage.xRange1)):
                self.slices.points1.append(f(self.odImage.xRange1[k], ch1[k])[0])
                self.slices.fit1.append(gaussian(resLSQ.x, [ch1[k], r[1][k]], 0)[0])
            self.slices.ch1 = ch1

        elif self.fitFunction == FIT_FUNCTIONS.index("Gaussian Mask Sigma"):
            # We mask EXCLUSION_RADIUS * initial sigma around the center of the fit
            EXCLUSION_RADIUS = self.WingRad

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D, Y2D = np.meshgrid(np.arange(dx), np.arange(dy))
                A = np.matrix(
                    np.column_stack((X2D.ravel(), Y2D.ravel(), np.ones(dx * dy)))
                )
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg = np.reshape(C * A.T, (dy, dx))
                return np.asarray(od - bg)

            # Gaussian fit with gradient TO BE UPDATED TO MASK FITTING
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg, 5, mode="constant")
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
                0,
                0,
            ]
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                2 * np.pi,
                40,
                40,
            ]
            pLower = [
                -np.inf,
                0.0,
                np.min(r[0]),
                0,
                np.min(r[1]),
                0,
                -2 * np.pi,
                -40,
                -40,
            ]
            p0 = checkGuess(p0, pUpper, pLower)

            pks = peak_local_max(blur, min_distance=20, exclude_border=2, num_peaks=3)
            guesses = [p0]
            for pk in pks:
                yc = pk[0]
                xc = pk[1]
                peak = data[yc, xc]
                offset = np.mean(data)
                if peak > 0:
                    (sigx, sigy) = 15, 15
                    guess = [
                        offset,
                        peak,
                        xmin + xc,
                        sigx,
                        ymin + yc,
                        sigy,
                        0,
                        0.0,
                        0.0,
                    ]
                    guesses.append(checkGuess(guess, pUpper, pLower))

            best_fit = None
            best_guess = np.inf

            # Initial Gaussian Fit
            for guess in guesses:
                try:
                    # print("Trying guess: {}".format(guess))
                    res = least_squares(
                        gaussian_mask_sigma,
                        guess,
                        args=(r, self.odImage.ODCorrected),
                        bounds=(pLower, pUpper),
                        kwargs={
                            "mask_radius_x": 1e-6,
                            "mask_radius_y": 1e-6,
                            "x0_mask": self.odImage.ODCorrected.shape[1] / 2,
                            "y0_mask": self.odImage.ODCorrected.shape[0] / 2,
                            "theta_mask": 0,
                        },
                    )
                    # print("Cost: {}".format(res.cost))
                    if not res.success:
                        print("Warning: fit did not converge.")
                    elif res.cost < best_guess:
                        best_guess = res.cost
                        best_fit = res
                except ValueError as e:
                    print("Error fitting image: {}".format(e))

            updated_guess = res.x
            best_fit = None
            best_guess = np.inf
            # print('UPDATED GUESS: {}'.format(updated_guess))
            try:
                # print("Trying guess: {}".format(guess))
                res = least_squares(
                    gaussian_mask_sigma,
                    guess,
                    args=(r, self.odImage.ODCorrected),
                    bounds=(pLower, pUpper),
                    kwargs={
                        "mask_radius_x": EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                        "mask_radius_y": EXCLUSION_RADIUS * np.abs(updated_guess[5]),
                        "x0_mask": np.abs(updated_guess[2]),
                        "y0_mask": np.abs(updated_guess[4]),
                        "theta_mask": np.abs(updated_guess[6]),
                    },
                )
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
                self.fittedImage = gaussian_mask_sigma(
                    resLSQ.x,
                    r,
                    0,
                    mask_radius_x=EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                    mask_radius_y=EXCLUSION_RADIUS * np.abs(updated_guess[5]),
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                    theta_mask=np.abs(updated_guess[6]),
                ).reshape(self.odImage.ODCorrected.shape)

                fittedImage_no_mask = gaussian_mask_sigma(
                    resLSQ.x,
                    r,
                    0,
                    mask_radius_x=1e-9,
                    mask_radius_y=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                    theta_mask=np.abs(updated_guess[6]),
                ).reshape(self.odImage.ODCorrected.shape)

            # Integration of number from computed column density

            # Compute the number by summing the pixels and multiplying by the pixel area
            # subtract the offset from the initial Gaussian fit, as well as the gradients
            # We make a background image where we just set the amplitude to zero, so that
            # we only have the offset and the gradients
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
            bg_res = [
                resLSQ.x[0],
                0,
                resLSQ.x[2],
                resLSQ.x[3],
                resLSQ.x[4],
                resLSQ.x[5],
                resLSQ.x[6],
                resLSQ.x[7],
                resLSQ.x[8],
            ]

            # fitted bg converted to density
            getsigma = SIGMA_0[self.odImage.species] * (
                2 if IMFIT_MODES[self.odImage.mode]["Image Path"] == "Vertical" else 1
            )

            fitted_bg = (
                gaussian_mask_sigma(
                    bg_res,
                    r,
                    0,
                    mask_radius_x=1e-9,
                    mask_radius_y=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                    theta_mask=np.abs(updated_guess[6]),
                ).reshape(self.odImage.ODCorrected.shape)
                / getsigma
            )

            raw_number = (self.odImage.n - fitted_bg).sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            error = (
                np.sqrt((self.odImage.nerr**2).sum())
                * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            self.fitData = np.append(self.fitData, number)
            self.fitData = np.append(self.fitData, error)
            self.fitData = np.append(self.fitData, EXCLUSION_RADIUS)

            ### Get radial average
            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]
            self.slices.fit0a = fittedImage_no_mask[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]
            self.slices.fit1a = fittedImage_no_mask[:, I0]

            print("Done with fit function!")

        elif self.fitFunction == FIT_FUNCTIONS.index("Gaussian Mask Sigma No Rot"):
            # We mask EXCLUSION_RADIUS * initial sigma around the center of the fit
            EXCLUSION_RADIUS = self.WingRad

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D, Y2D = np.meshgrid(np.arange(dx), np.arange(dy))
                A = np.matrix(
                    np.column_stack((X2D.ravel(), Y2D.ravel(), np.ones(dx * dy)))
                )
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg = np.reshape(C * A.T, (dy, dx))
                return np.asarray(od - bg)

            # Gaussian fit with gradient TO BE UPDATED TO MASK FITTING
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg, 5, mode="constant")
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
                0,
            ]
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                40,
                40,
            ]
            pLower = [
                -np.inf,
                0.0,
                np.min(r[0]),
                0,
                np.min(r[1]),
                0,
                -40,
                -40,
            ]
            p0 = checkGuess(p0, pUpper, pLower)

            pks = peak_local_max(blur, min_distance=20, exclude_border=2, num_peaks=3)
            guesses = [p0]
            for pk in pks:
                yc = pk[0]
                xc = pk[1]
                peak = data[yc, xc]
                offset = np.mean(data)
                if peak > 0:
                    (sigx, sigy) = 15, 15
                    guess = [offset, peak, xmin + xc, sigx, ymin + yc, sigy, 0.0, 0.0]
                    guesses.append(checkGuess(guess, pUpper, pLower))

            best_fit = None
            best_guess = np.inf

            # Initial Gaussian Fit
            for guess in guesses:
                try:
                    # print("Trying guess: {}".format(guess))
                    res = least_squares(
                        gaussian_mask_sigma_no_rot,
                        guess,
                        args=(r, self.odImage.ODCorrected),
                        bounds=(pLower, pUpper),
                        kwargs={
                            "mask_radius_x": 1e-6,
                            "mask_radius_y": 1e-6,
                            "x0_mask": self.odImage.ODCorrected.shape[1] / 2,
                            "y0_mask": self.odImage.ODCorrected.shape[0] / 2,
                        },
                    )
                    # print("Cost: {}".format(res.cost))
                    if not res.success:
                        print("Warning: fit did not converge.")
                    elif res.cost < best_guess:
                        best_guess = res.cost
                        best_fit = res
                except ValueError as e:
                    print("Error fitting image: {}".format(e))

            updated_guess = res.x
            best_fit = None
            best_guess = np.inf
            # print('UPDATED GUESS: {}'.format(updated_guess))
            try:
                # print("Trying guess: {}".format(guess))
                res = least_squares(
                    gaussian_mask_sigma_no_rot,
                    guess,
                    args=(r, self.odImage.ODCorrected),
                    bounds=(pLower, pUpper),
                    kwargs={
                        "mask_radius_x": EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                        "mask_radius_y": EXCLUSION_RADIUS * np.abs(updated_guess[5]),
                        "x0_mask": np.abs(updated_guess[2]),
                        "y0_mask": np.abs(updated_guess[4]),
                    },
                )
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
                self.fittedImage = gaussian_mask_sigma_no_rot(
                    resLSQ.x,
                    r,
                    0,
                    mask_radius_x=EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                    mask_radius_y=EXCLUSION_RADIUS * np.abs(updated_guess[5]),
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                ).reshape(self.odImage.ODCorrected.shape)
                # For plotting:
                fittedImage_no_mask = gaussian_mask_sigma_no_rot(
                    resLSQ.x,
                    r,
                    0,
                    mask_radius_x=1e-9,
                    mask_radius_y=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                ).reshape(self.odImage.ODCorrected.shape)

            # Integration of number from computed column density

            # Compute the number by summing the pixels and multiplying by the pixel area
            # subtract the offset from the initial Gaussian fit, as well as the gradients
            # We make a background image where we just set the amplitude to zero, so that
            # we only have the offset and the gradients
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
            bg_res = [
                resLSQ.x[0],
                0,
                resLSQ.x[2],
                resLSQ.x[3],
                resLSQ.x[4],
                resLSQ.x[5],
                resLSQ.x[6],
                resLSQ.x[7],
            ]

            # fitted bg converted to density
            getsigma = SIGMA_0[self.odImage.species] * (
                2 if IMFIT_MODES[self.odImage.mode]["Image Path"] == "Vertical" else 1
            )

            fitted_bg = (
                gaussian_mask_sigma_no_rot(
                    bg_res,
                    r,
                    0,
                    mask_radius_x=1e-9,
                    mask_radius_y=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                ).reshape(self.odImage.ODCorrected.shape)
                / getsigma
            )

            raw_number = (self.odImage.n - fitted_bg).sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            error = (
                np.sqrt((self.odImage.nerr**2).sum())
                * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            self.fitData = np.append(self.fitData, number)
            self.fitData = np.append(self.fitData, error)
            self.fitData = np.append(self.fitData, EXCLUSION_RADIUS)

            ### Get radial average
            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]
            self.slices.fit0a = fittedImage_no_mask[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]
            self.slices.fit1a = fittedImage_no_mask[:, I0]

            print("Done with fit function!")

        elif self.fitFunction == FIT_FUNCTIONS.index("Twisted Gaussian"):

            # Gaussian fit, rotated by VERT_TRAP_ANGLE

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]
            pUpper = [np.inf, 50.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0, pUpper, pLower)

            angle = IMFIT_MODES[self.mode]["Fit angle"]
            resLSQ = least_squares(
                gaussianNoRotTwist,
                p0,
                args=(r, self.odImage.ODCorrected, angle),
                bounds=(pLower, pUpper),
            )

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussianNoRotTwist(resLSQ.x, r, 0, angle).reshape(
                self.odImage.ODCorrected.shape
            )

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ####### Calculate slices through fit #####

            f = RectBivariateSpline(
                self.odImage.xRange1,
                self.odImage.xRange0,
                self.odImage.ODCorrected,
            )

            m0 = np.tan(np.pi / 2.0 - angle * np.pi / 180.0)
            m1 = -np.tan(angle * np.pi / 180.0)

            if abs(m0) > abs(m1):
                # Ensure the that lower slope is always along x
                m0, m1 = m1, m0

            b0 = -m0 * resLSQ.x[2] + resLSQ.x[4]
            ch0 = np.asarray(self.odImage.xRange0) * m0 + b0
            b1 = -m1 * resLSQ.x[2] + resLSQ.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1) / m1

            for k in range(len(self.odImage.xRange0)):
                self.slices.points0.append(f(ch0[k], self.odImage.xRange0[k])[0])
                self.slices.fit0.append(
                    gaussianNoRotTwist(resLSQ.x, [r[0][k], ch0[k]], 0.0, angle)[0]
                )
            self.slices.ch0 = ch0

            for k in range(len(self.odImage.xRange1)):
                self.slices.points1.append(f(self.odImage.xRange1[k], ch1[k])[0])
                self.slices.fit1.append(
                    gaussianNoRotTwist(resLSQ.x, [ch1[k], r[1][k]], 0, angle)[0]
                )
            self.slices.ch1 = ch1

        elif self.fitFunction == FIT_FUNCTIONS.index("Gaussian"):

            # Gaussian fit without rotation

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]
            pUpper = [np.inf, 50.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                gaussianNoRot,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = gaussianNoRot(resLSQ.x, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

        elif self.fitFunction == FIT_FUNCTIONS.index("Gaussian w/ Gradient"):

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D, Y2D = np.meshgrid(np.arange(dx), np.arange(dy))
                A = np.matrix(
                    np.column_stack((X2D.ravel(), Y2D.ravel(), np.ones(dx * dy)))
                )
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg = np.reshape(C * A.T, (dy, dx))
                return np.asarray(od - bg)

            # Gaussian fit without rotation
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg, 5, mode="constant")

            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
                0,
            ]
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                40,
                40,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -40, -40]
            p0 = checkGuess(p0, pUpper, pLower)

            pks = peak_local_max(blur, min_distance=20, exclude_border=2, num_peaks=3)
            guesses = [p0]
            for pk in pks:
                yc = pk[0]
                xc = pk[1]
                peak = data[yc, xc]
                offset = np.mean(data)
                if peak > 0:
                    (sigx, sigy) = 15, 15
                    guess = [offset, peak, xmin + xc, sigx, ymin + yc, sigy, 0.0, 0.0]
                    guesses.append(checkGuess(guess, pUpper, pLower))

            best_fit = None
            best_guess = np.inf
            for guess in guesses:
                try:
                    # print("Trying guess: {}".format(guess))
                    res = least_squares(
                        gaussianNoRotGradient,
                        guess,
                        args=(r, self.odImage.ODCorrected),
                        bounds=(pLower, pUpper),
                    )
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
                self.fittedImage = gaussianNoRotGradient(resLSQ.x, r, 0).reshape(
                    self.odImage.ODCorrected.shape
                )

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

            print("Done with fit function!")

        elif self.fitFunction == FIT_FUNCTIONS.index("Gaussian Fixed"):

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D, Y2D = np.meshgrid(np.arange(dx), np.arange(dy))
                A = np.matrix(
                    np.column_stack((X2D.ravel(), Y2D.ravel(), np.ones(dx * dy)))
                )
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg = np.reshape(C * A.T, (dy, dx))
                return np.asarray(od - bg)

            # Gaussian fit without rotation
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])
            xmax = np.max(r[0])
            ymax = np.max(r[1])
            xc = 0.5 * (xmin + xmax)
            yc = 0.5 * (ymin + ymax)

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg, 5, mode="constant")

            p0 = [0, M, xc, 15, yc, 3, 0, 0]
            pUpper = [np.inf, 50.0, xc + 2, 40, yc + 2, 6, 40, 40]
            pLower = [-np.inf, 0.0, xc - 2, 4, yc - 2, 0.01, -40, -40]
            p0 = checkGuess(p0, pUpper, pLower)

            guesses = [p0]

            best_fit = None
            best_guess = np.inf
            for guess in guesses:
                try:
                    print("Trying guess: {}".format(guess))
                    res = least_squares(
                        gaussianNoRotGradient,
                        guess,
                        args=(r, self.odImage.ODCorrected),
                        bounds=(pLower, pUpper),
                    )
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
                self.fittedImage = gaussianNoRotGradient(resLSQ.x, r, 0).reshape(
                    self.odImage.ODCorrected.shape
                )

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

            print("Done with fit function!")

        elif self.fitFunction == FIT_FUNCTIONS.index("Bigaussian"):
            # Bi Gaussian fit
            ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]
            p0 = [
                0,
                M / 2.0,
                20.0,
                20.0,
                M / 2.0,
                8.0,
                8.0,
                self.odImage.xRange0[I1],
                self.odImage.xRange1[I0],
            ]
            pUpper = [
                np.inf,
                8.0,
                len(r[0]),
                len(r[1]),
                8.0,
                len(r[0]),
                len(r[1]),
                np.max(r[0]),
                np.max(r[1]),
            ]
            pLower = [-np.inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                doubleGaussian,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )
            # self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = doubleGaussian(resLSQ.x, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            ### Make sure the BEC comes first in the list

            if self.fitData[2] + self.fitData[3] > self.fitData[5] + self.fitData[6]:
                t = copy.deepcopy(self.fitData[4:7])

                self.fitData[4:7] = self.fitData[1:4]
                self.fitData[1:4] = t
            else:
                pass

            ### Calculate radial average
            x_coerced = min(
                max(self.fitData[7], self.odImage.xRange0.start),
                self.odImage.xRange0.stop - 1,
            )
            y_coerced = min(
                max(self.fitData[8], self.odImage.xRange1.start),
                self.odImage.xRange1.stop - 1,
            )

            I0 = self.odImage.xRange0.index(int(x_coerced))
            I1 = self.odImage.xRange1.index(int(y_coerced))

            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, [I0, I1])
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, [I0, I1])
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

        elif self.fitFunction == FIT_FUNCTIONS.index("Thomas-Fermi"):
            # Thomas-Fermi fit
            ### Parameters: [offset, ampTF, x0, rx, y0, ry, ampGauss, wx, wy]
            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            thermal_sizes = [5, 10, 20, 40, 80]
            resLSQ = None
            for ts in thermal_sizes:
                p0 = [
                    0,
                    M,
                    self.odImage.xRange0[I1],
                    10,
                    self.odImage.xRange1[I0],
                    10,
                    0.1,
                    ts,
                    ts,
                ]
                pUpper = [
                    np.inf,
                    50.0,
                    np.max(r[0]),
                    len(r[0]),
                    np.max(r[1]),
                    len(r[1]),
                    50.0,
                    len(r[0]),
                    len(r[1]),
                ]
                pLower = [-np.inf, 0, np.min(r[0]), 0, np.min(r[1]), 0, 0, 0, 0]
                p0 = checkGuess(p0, pUpper, pLower)
                res = least_squares(
                    thomasFermi,
                    p0,
                    args=(r, self.odImage.ODCorrected),
                    bounds=(pLower, pUpper),
                )
                if resLSQ is None or res.cost < resLSQ.cost:
                    resLSQ = res
            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = thomasFermi(resLSQ.x, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )
            # T-F: [offset, ampTF, x0, rx, y0, ry, ampGauss, wx, wy]
            # Gaussian: [offset, amplitude, x0, wx, y0, w]
            self.fittedImageGauss = gaussianNoRot(
                [
                    resLSQ.x[0],
                    resLSQ.x[6],
                    resLSQ.x[2],
                    resLSQ.x[7],
                    resLSQ.x[4],
                    resLSQ.x[8],
                ],
                r,
                0,
            ).reshape(self.odImage.ODCorrected.shape)

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit

            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]
            self.slices.fit0Gauss = self.fittedImageGauss[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]
            self.slices.fit1Gauss = self.fittedImageGauss[:, I0]

        elif self.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac"):

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            # Fermi--Dirac fit
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]
            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
            ]  # An initial q of 0 corresponds to T/TF=0.56
            pUpper = [
                np.inf,
                900.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                np.inf,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                fermiDirac,
                p0,
                args=(r, self.odImage.ODCorrected),
                kwargs={"mask_above": MAX_OD_FIT},
                bounds=(pLower, pUpper),
            )
            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = fermiDirac(resLSQ.x, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            # Gaussian fit
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0]

            pUpper = [
                np.inf,
                100.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                0.00001,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -0.00001]

            resLSQ = least_squares(
                gaussian,
                p0,
                args=(r, self.odImage.ODCorrected),
                kwargs={"mask_above": MAX_OD_FIT},
                bounds=(pLower, pUpper),
            )
            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussian(self.fitDataGauss, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            #######################################################################################################################

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            # azAverage -- I'm not sure what the correct index should be for azAverage...
            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = azimuthalAverage(
                self.fittedImageGauss, center
            )

            ### Calculate slices through fit (No rotation)
            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

        elif self.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac fixed betamu"):
            ### Parameters: [offset, amplitude, x0, wx, y0, wy]

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            # INITIAL Gaussian fit with gradient
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
                0,
                0,
            ]
            # we barely allow an angle
            pUpper = [
                np.inf,
                100.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                0.00001,
                np.inf,
                np.inf,
            ]
            pLower = [
                -np.inf,
                0.0,
                np.min(r[0]),
                0,
                np.min(r[1]),
                0,
                -0.00001,
                -np.inf,
                -np.inf,
            ]

            resLSQ = least_squares(
                gaussianGradient,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )

            # Integration of number from computed column density. This is to give
            # a relation between T and betamu since N = -(k_B T / (hbar * omega_bar))^3 Li_3(-exp(betamu))
            # we take it as a constraint for our fit

            # Compute the number by summing the pixels and multiplying by the pixel area
            # subtract the offset from the initial Gaussian fit, as well as the gradients
            # We make a background image where we just set the amplitude to zero, so that
            # we only have the offset and the gradients
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
            bg_res = [
                resLSQ.x[0],
                0,
                resLSQ.x[2],
                resLSQ.x[3],
                resLSQ.x[4],
                resLSQ.x[5],
                resLSQ.x[6],
                resLSQ.x[7],
                resLSQ.x[8],
            ]

            # fitted bg converted to density
            getsigma = SIGMA_0[self.odImage.species] * (
                2 if IMFIT_MODES[self.odImage.mode]["Image Path"] == "Vertical" else 1
            )

            fitted_bg = (
                gaussianGradient(
                    bg_res,
                    r,
                    0,
                ).reshape(self.odImage.ODCorrected.shape)
                / getsigma
            )

            raw_number = (self.odImage.n - fitted_bg).sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )
            print("Computed number from image: {:.2e}".format(number))
            
            # multiply the number with the fudge factor
            number *= self.Nscaler

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            # Fermi--Dirac fit with fixed betamu
            ### Parameters: [offset, amplitude, x0, wx, y0, wy]
            p0 = [
                resLSQ.x[0],
                resLSQ.x[1],
                resLSQ.x[2],
                resLSQ.x[3],
                resLSQ.x[4],
                resLSQ.x[5],
            ]  # An initial q of 0 corresponds to T/TF=0.56
            pUpper = [
                np.inf,
                np.inf,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0, pUpper, pLower)

            TOF0 = self.TOF * 1e-3  # Convert TOF from ms to s


            # Other input constants
            hbar =  NAT_CONSTANTS["hbar"]  # J s
            kB = NAT_CONSTANTS["kB"]  # J/K
            omega_x = 2 * np.pi * self.fx
            omega_y = 2 * np.pi * self.fy
            omega_z = 2 * np.pi * self.fz
            amu2kg = NAT_CONSTANTS["amu2kg"]  # kg
            # mass = 40 * amu2kg

            if self.mbemu == "K":
                mass0 = 40
                mass = mass0 * amu2kg  # kg
            elif self.mbemu == "KRb":
                mass0 = 127
                mass = mass0 * amu2kg  # kg
            else:
                raise ValueError("Unknown species: {}".format(self.odImage.species))
            
            pxsz_um = self.config["Pixel Size"] * self.odImage.data.bin * 1e-6
            print("Pixel size: {:f} um".format(pxsz_um * 1e6))

            resLSQ = least_squares(
                fermiDirac_fixed_bemu,
                p0,
                args=(r, self.odImage.ODCorrected),
                kwargs={
                    "N0": number,
                    "TOF": TOF0,
                    "omega_x": omega_x,
                    "omega_y": omega_y,
                    "omega_z": omega_z,
                    "pxsz_um": pxsz_um,
                    "mass": mass,
                },  # Constrain THE NUMBER
                bounds=(pLower, pUpper),
                xtol=3e-16,  # Change tolerance for the fit; needed for convergence
                ftol=3e-16,
                diff_step=[1e-12, 1e-12, 0.01, 0.01, 0.01, 0.01],
            )

            # From fit get betamu (needs to be cleaned once we use the inputs)

            omega_bar = (omega_x * omega_y * omega_z) ** (1 / 3)
            Tx = (
                mass
                * omega_x**2
                * (resLSQ.x[3] * pxsz_um) ** 2
                / (1 + omega_x**2 * TOF0**2)
                / kB
            )
            Ty = (
                mass
                * omega_y**2
                * (resLSQ.x[5] * pxsz_um) ** 2
                / (1 + omega_y**2 * TOF0**2)
                / kB
            )
            T_avg = (Tx**2 * Ty) ** (1 / 3)
            print("Fitted T_avg: {:f} nK".format(T_avg * 1e9))

            # Find where betamu for the given T gives us the right number of particles
            betamu_range = np.linspace(-10, 20, 5000)
            N_checker_3D = (
                (kB * T_avg / (hbar * omega_bar)) ** 3
            ) * polylog_lib.fermi_poly3(betamu_range)

            N_diff_3D = np.abs(N_checker_3D - number)
            ind_N_3D = np.where(N_diff_3D == np.min(N_diff_3D))[0]
            betamu_3D = betamu_range[ind_N_3D[0]]

            print("Computed betamu: {:.2f}".format(betamu_3D))

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fitData = np.append(self.fitData, betamu_3D)
            self.fitData = np.append(self.fitData, self.TOF)
            self.fitData = np.append(self.fitData, self.fx)
            self.fitData = np.append(self.fitData, self.fy)
            self.fitData = np.append(self.fitData, self.fz)
            self.fitData = np.append(self.fitData, number)
            self.fitData = np.append(self.fitData, self.Nscaler)
            self.fitData = np.append(self.fitData, mass0)

            self.fittedImage = fermiDirac_fixed_bemu(
                resLSQ.x,
                r,
                0,
                N0=number,
                TOF=TOF0,
                omega_x=omega_x,
                omega_y=omega_y,
                omega_z=omega_z,
                pxsz_um = pxsz_um,
                mass = mass,
            ).reshape(self.odImage.ODCorrected.shape)

            # Gaussian fit
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20, 0]

            pUpper = [
                np.inf,
                100.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                0.00001,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -0.00001]

            resLSQ = least_squares(
                gaussian,
                p0,
                args=(r, self.odImage.ODCorrected),
                kwargs={"mask_above": MAX_OD_FIT},
                bounds=(pLower, pUpper),
            )
            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussian(self.fitDataGauss, r, 0).reshape(
                self.odImage.ODCorrected.shape
            )

            #######################################################################################################################

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            # azAverage -- I'm not sure what the correct index should be for azAverage...
            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = azimuthalAverage(
                self.fittedImageGauss, center
            )

            ### Calculate slices through fit (No rotation)
            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = self.fittedImage[I1, :]

            self.slices.points1 = self.odImage.ODCorrected[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = self.fittedImage[:, I0]

        elif self.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac 2D"):

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            angle = IMFIT_MODES[self.mode]["Fit angle"]

            # Gaussian fit
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]

            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]

            resLSQ = least_squares(
                gaussianNoRotTwist,
                p0,
                args=(r, self.odImage.ODCorrected, angle),
                bounds=(pLower, pUpper),
            )
            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussianNoRotTwist(
                self.fitDataGauss, r, 0, angle
            ).reshape(self.odImage.ODCorrected.shape)

            print("Gaussian fit parameters: {}".format(self.fitDataGauss))

            # Fermi--Dirac fit
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]
            p0 = [
                self.fitDataGauss[0],
                self.fitDataGauss[1],
                self.fitDataGauss[2],
                self.fitDataGauss[3],
                self.fitDataGauss[4],
                self.fitDataGauss[5],
                0,
            ]  # An initial q of 0 corresponds to T/TF=0.56
            pUpper = [
                np.inf,
                900.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                np.inf,
            ]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                fermiDirac2D,
                p0,
                args=(r, self.odImage.ODCorrected, angle),
                bounds=(pLower, pUpper),
            )
            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = fermiDirac2D(self.fitData, r, 0, angle).reshape(
                self.odImage.ODCorrected.shape
            )

            print("Fermi-Dirac fit parameters: {}".format(self.fitData))

            ### Get radial average

            I0 = self.odImage.xRange0.index(int(self.fitData[2]))
            I1 = self.odImage.xRange1.index(int(self.fitData[4]))

            center = [I0, I1]
            self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            self.slices.radSliceFitGauss = azimuthalAverage(
                self.fittedImageGauss, center
            )

            ####### Calculate slices through fit #####

            f = RectBivariateSpline(
                self.odImage.xRange1,
                self.odImage.xRange0,
                self.odImage.ODCorrected,
            )

            m0 = np.tan(np.pi / 2.0 - angle * np.pi / 180.0)
            m1 = -np.tan(angle * np.pi / 180.0)

            if abs(m0) > abs(m1):
                # Ensure the that lower slope is always along x
                m0, m1 = m1, m0

            b0 = -m0 * resLSQ.x[2] + resLSQ.x[4]
            ch0 = np.asarray(self.odImage.xRange0) * m0 + b0
            b1 = -m1 * resLSQ.x[2] + resLSQ.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1) / m1

            ### Calculate slices through fit
            for k in range(len(self.odImage.xRange0)):
                self.slices.points0.append(f(ch0[k], self.odImage.xRange0[k])[0])
                self.slices.fit0.append(
                    # gaussianNoRotTwist(resLSQ.x, [r[0][k], ch0[k]], 0.0, angle)[0]
                    fermiDirac2D(resLSQ.x, np.array([r[0][k], ch0[k]]), 0, angle)[0]
                )
            self.slices.ch0 = ch0

            for k in range(len(self.odImage.xRange1)):
                self.slices.points1.append(f(self.odImage.xRange1[k], ch1[k])[0])
                self.slices.fit1.append(
                    # gaussianNoRotTwist(resLSQ.x, [ch1[k], r[1][k]], 0, angle)[0]
                    fermiDirac2D(resLSQ.x, np.array([ch1[k], r[1][k]]), 0, angle)[0]
                )
            self.slices.ch1 = ch1

        elif self.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac 2D Int"):

            r = [None, None]
            r[0] = self.odImage.xRange0  # x
            r[1] = self.odImage.xRange1  # y

            # Fit integrated OD
            od_int = np.sum(self.odImage.ODCorrected, axis=0)
            od_int *= self.config["Pixel Size"] * self.odImage.data.bin

            # subtract gradient using linear fit
            A = np.vstack([np.arange(len(od_int)), np.ones(len(od_int))]).T
            m, c = np.linalg.lstsq(A, od_int, rcond=None)[0]
            od_int_no_bg = od_int - (m * np.arange(len(od_int)) + c)

            # 2D Gaussian fit
            # Parameters: [offset, amplitude, x0, wx, y0, wy]
            # Only to get the center point for the x slice
            p0 = [0, M, self.odImage.xRange0[I1], 20, self.odImage.xRange1[I0], 20]
            pUpper = [np.inf, 50.0, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, np.min(r[0]), 0, np.min(r[1]), 0]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                gaussianNoRot,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )

            y_for_x_slice = resLSQ.x[4] - r[1][0]
            x_for_y_slice = resLSQ.x[2] - r[0][0]

            # 1D Gaussian fit to integrated OD
            # Parameters: [offset, amplitude, x0, wx]
            est_bg = np.quantile(od_int, 0.1)
            p0 = [
                est_bg,
                np.max(od_int_no_bg),
                np.argmax(od_int_no_bg) + r[0][0],
                np.std(od_int_no_bg),
                m,
            ]
            pUpper = [
                np.inf,
                np.max(od_int) - np.min(od_int),
                r[0][-1],
                len(r[0]),
                np.inf,
            ]
            pLower = [-np.inf, 0.0, r[0][0], 0.0, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                gaussian1D,
                p0,
                args=(np.array(r[0]), od_int),
                bounds=(pLower, pUpper),
            )

            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussian1D(resLSQ.x, np.array(r[0]), 0)

            self.slices.points0 = od_int
            self.slices.points1 = self.odImage.ODCorrected[
                round(y_for_x_slice), :  # bug fixed 2025/05/28; vertical slice
            ]
            self.slices.ch0 = np.ones(len(r[0])) * y_for_x_slice + r[1][0]
            self.slices.ch1 = None
            self.slices.fit0 = self.fittedImageGauss

            # Fermi--Dirac fit
            # Parameters: [offset, amplitude, x0, sigma, q, gradient]
            p0 = [
                resLSQ.x[0],
                resLSQ.x[1],
                resLSQ.x[2],
                resLSQ.x[3],
                0,
                resLSQ.x[4],
            ]
            pUpper = [
                np.inf,
                np.inf,
                r[0][-1],
                len(r[0]),
                np.inf,
                np.inf,
            ]
            pLower = [-np.inf, 0.0, r[0][0], 0.0, -np.inf, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ = least_squares(
                fermiDirac2Dint,
                p0,
                args=(np.array(r[0]), od_int),
                bounds=(pLower, pUpper),
            )

            self.fitDataConf = confidenceIntervals(resLSQ)
            self.fitData = resLSQ.x
            self.fittedImage = fermiDirac2Dint(resLSQ.x, np.array(r[0]), 0)

            self.slices.fit1 = self.fittedImage

        elif self.fitFunction == FIT_FUNCTIONS.index("Twisted Fermi-Dirac 2D Int"):

            r = [None, None]
            r[0] = self.odImage.xRange0  # x
            r[1] = self.odImage.xRange1  # y

            angle = IMFIT_MODES[self.mode]["Fit angle"]
            # Rotate the image
            rot_img = ndimage.rotate(self.odImage.ODCorrected,
                                      -IMFIT_MODES[self.mode]["Fit angle"],
                                      reshape=False,)
                                   
            r_rot = [None, None]
            r_rot[0] = np.arange(0, np.shape(rot_img)[1])  # x
            r_rot[1] = np.arange(0, np.shape(rot_img)[0])  # y

            # Integrate the OD along x and y to get 1D profiles
            od_int_x = np.sum(rot_img, axis = 1)
            od_int_y = np.sum(rot_img, axis = 0)

            # Fit the 1D profiles with Fermi-Dirac 2D Int

            # Fit integrated OD
            od_int_x *= self.config["Pixel Size"] * self.odImage.data.bin
            od_int_y *= self.config["Pixel Size"] * self.odImage.data.bin

            # subtract gradient using linear fit
            A = np.vstack([np.arange(len(od_int_x)), np.ones(len(od_int_x))]).T
            m, c = np.linalg.lstsq(A, od_int_x, rcond=None)[0]
            od_int_x_no_bg = od_int_x - (m * np.arange(len(od_int_x)) + c)
            
            A = np.vstack([np.arange(len(od_int_y)), np.ones(len(od_int_y))]).T
            m, c = np.linalg.lstsq(A, od_int_y, rcond=None)[0]
            od_int_y_no_bg = od_int_y - (m * np.arange(len(od_int_y)) + c)

            # 2D Gaussian fit
            ### Parameters: [offset, amplitude, x0, wx, y0, wy]
            # Only to get the center point for the x slice
            p0 = [0, M, r[0][0] + np.shape(self.odImage.ODCorrected)[1] // 2, 20,
                  r[1][0] + np.shape(self.odImage.ODCorrected)[0] // 2, 20]
            pUpper = [np.inf, 50, np.max(r[0]), len(r[0]), np.max(r[1]), len(r[1])]
            pLower = [-np.inf, 0.0, 0, 0, 0, 0]

            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ_G_norot = least_squares(
                gaussianNoRotTwist,
                p0,
                args=(r, self.odImage.ODCorrected, angle),
                bounds=(pLower, pUpper),
            )

            self.fitDataConfGauss = confidenceIntervals(resLSQ_G_norot)
            self.fitDataGauss = resLSQ_G_norot.x
            self.fittedImageGauss = gaussianNoRotTwist(resLSQ_G_norot.x, r, 0, angle).reshape(self.odImage.ODCorrected.shape)

            ####### Calculate slices through fit #####
            m0 = np.tan(np.pi / 2.0 - angle * np.pi / 180.0)
            m1 = -np.tan(angle * np.pi / 180.0)

            if abs(m0) > abs(m1):
                # Ensure that lower slope is always along x
                m0, m1 = m1, m0

            b0 = -m0 * resLSQ_G_norot.x[2] + resLSQ_G_norot.x[4]
            ch0 = np.asarray(self.odImage.xRange0) * m0 + b0
            b1 = -m1 * resLSQ_G_norot.x[2] + resLSQ_G_norot.x[4]
            ch1 = (np.asarray(self.odImage.xRange1) - b1) / m1

            # 1D Gaussian fit along x to integrated OD along y
            # Parameters: [offset, amplitude, x0, wx]
            est_bg = np.quantile(od_int_y_no_bg, 0.1)
            p0 = [
                est_bg,
                np.max(od_int_y_no_bg),
                np.argmax(od_int_y_no_bg) + r_rot[0][0],
                np.std(od_int_y_no_bg),
                m,
            ]
            pUpper = [
                np.inf,
                np.max(od_int_y_no_bg) - np.min(od_int_y_no_bg),
                r_rot[0][-1],
                len(r_rot[0]),
                np.inf,
            ]
            pLower = [-np.inf, 0.0, r_rot[0][0], 0.0, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ_Gx = least_squares(
                gaussian1D,
                p0,
                args=(np.array(r_rot[0]), od_int_y_no_bg),
                bounds=(pLower, pUpper),
            )

            ### Parameters: [offset, amplitude, x0, wx, y0, wy]
            self.slices.points0 = od_int_y_no_bg
            self.slices.points1 = od_int_x_no_bg
            self.slices.ch0 = ch0
            self.slices.ch1 = ch1

            # Fermi--Dirac fit
            # Parameters: [offset, amplitude, x0, sigma, q, gradient]
            p0 = [
                resLSQ_Gx.x[0],
                resLSQ_Gx.x[1],
                resLSQ_Gx.x[2],
                resLSQ_Gx.x[3],
                0,
                resLSQ_Gx.x[4],
            ]
            pUpper = [
                np.inf,
                np.inf,
                r_rot[1][-1],
                len(r_rot[0]),
                np.inf,
                np.inf,
            ]
            pLower = [-np.inf, 0.0, r_rot[1][0], 0.0, -np.inf, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ_inty = least_squares(
                fermiDirac2Dint,
                p0,
                args=(np.array(r_rot[0]), od_int_y_no_bg),
                bounds=(pLower, pUpper),
            )

            self.fittedImage = fermiDirac2Dint(resLSQ_inty.x, np.array(r_rot[0]), 0)
            self.slices.fit0 = self.fittedImage

            # 1D Gaussian fit along y to integrated OD along x
            resLSQ_Gy = least_squares(
                gaussian1D,
                p0,
                args=(np.array(r_rot[1]), od_int_x_no_bg),
                bounds=(pLower, pUpper),
            )

            # Fermi--Dirac fit
            # Parameters: [offset, amplitude, x0, sigma, q, gradient]
            p0 = [
                resLSQ_Gy.x[0],
                resLSQ_Gy.x[1],
                resLSQ_Gy.x[2],
                resLSQ_Gy.x[3],
                0,
                resLSQ_Gy.x[4],
            ]
            pUpper = [
                np.inf,
                np.inf,
                r_rot[1][-1],
                len(r_rot[0]),
                np.inf,
                np.inf,
            ]
            pLower = [-np.inf, 0.0, r_rot[1][0], 0.0, -np.inf, -np.inf]
            p0 = checkGuess(p0, pUpper, pLower)

            resLSQ_intx = least_squares(
                fermiDirac2Dint,
                p0,
                args=(np.array(r_rot[1]), od_int_x_no_bg),
                bounds=(pLower, pUpper),
            )


            fitData = resLSQ_inty.x
            fitData = np.append(fitData, resLSQ_intx.x)

            # Calculate average number density in border and subtract from rest of image
            border = int(max(min(self.odImage.n.shape) / 10, 5))
            border_mask = np.ones(self.odImage.n.shape)
            border_mask[border:-border, border:-border] = 0
            offset = np.sum(self.odImage.n * border_mask) / np.sum(border_mask)
            self.odImage.n -= offset

            interior = self.odImage.n[border:-border, border:-border]
            interior_err = self.odImage.nerr[border:-border, border:-border]

            # Compute the number by summing the pixels and multiplying by the pixel area
            raw_number = interior.sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            error = (
                np.sqrt((interior_err**2).sum())
                * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )
            
            fitData = np.append(fitData, number)

            conf_int = confidenceIntervals(resLSQ_inty)
            conf_int = np.append(conf_int, confidenceIntervals(resLSQ_intx))

            self.fitDataConf = conf_int
            self.fitData = fitData
            self.fittedImage = fermiDirac2Dint(resLSQ_intx.x, np.array(r_rot[1]), 0)
            self.slices.fit1 = self.fittedImage


        elif self.fitFunction == FIT_FUNCTIONS.index("Gauss (Mask) Int"):
            # Gaussian fit with gradient in 1D with masked center
            ### Parameters: [offset, amplitude, x0, wx, dODdx]

            r = [None, None]
            r[0] = self.odImage.xRange0  # x
            r[1] = self.odImage.xRange1  # y

            # Fit integrated OD
            od_int = np.sum(self.odImage.ODCorrected, axis=0)
            od_int *= self.config["Pixel Size"] * self.odImage.data.bin

             # We mask EXCLUSION_RADIUS * initial sigma around the center of the fit
            EXCLUSION_RADIUS = self.WingRad

            def subtract_gradient(od):
                (dy, dx) = np.shape(od)
                X2D, Y2D = np.meshgrid(np.arange(dx), np.arange(dy))
                A = np.matrix(
                    np.column_stack((X2D.ravel(), Y2D.ravel(), np.ones(dx * dy)))
                )
                B = od.flatten()
                C = np.dot((A.T * A).I * A.T, B).flatten()
                bg = np.reshape(C * A.T, (dy, dx))
                return np.asarray(od - bg)

            # First perform a 2D Gaussian fit to get initial parameters
            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1
            xmin = np.min(r[0])
            ymin = np.min(r[1])

            data = self.odImage.ODCorrected
            od_no_bg = subtract_gradient(data)
            blur = ndimage.gaussian_filter(od_no_bg, 5, mode="constant")
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
            p0 = [
                0,
                M,
                self.odImage.xRange0[I1],
                20,
                self.odImage.xRange1[I0],
                20,
                0,
                0,
            ]
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                40,
                40,
            ]
            pLower = [
                -np.inf,
                0.0,
                np.min(r[0]),
                0,
                np.min(r[1]),
                0,
                -40,
                -40,
            ]
            p0 = checkGuess(p0, pUpper, pLower)

            pks = peak_local_max(blur, min_distance=20, exclude_border=2, num_peaks=3)
            guesses = [p0]
            for pk in pks:
                yc = pk[0]
                xc = pk[1]
                peak = data[yc, xc]
                offset = np.mean(data)
                if peak > 0:
                    (sigx, sigy) = 15, 15
                    guess = [offset, peak, xmin + xc, sigx, ymin + yc, sigy, 0.0, 0.0]
                    guesses.append(checkGuess(guess, pUpper, pLower))

            best_fit = None
            best_guess = np.inf

            # Initial Gaussian Fit
            for guess in guesses:
                try:
                    # print("Trying guess: {}".format(guess))
                    res = least_squares(
                        gaussian_mask_sigma_no_rot,
                        guess,
                        args=(r, self.odImage.ODCorrected),
                        bounds=(pLower, pUpper),
                        kwargs={
                            "mask_radius_x": 1e-6,
                            "mask_radius_y": 1e-6,
                            "x0_mask": self.odImage.ODCorrected.shape[1] / 2,
                            "y0_mask": self.odImage.ODCorrected.shape[0] / 2,
                        },
                    )
                    # print("Cost: {}".format(res.cost))
                    if not res.success:
                        print("Warning: fit did not converge.")
                    elif res.cost < best_guess:
                        best_guess = res.cost
                        best_fit = res
                except ValueError as e:
                    print("Error fitting image: {}".format(e))
            res2D = res.x

            updated_guess = res.x[[0,1,2,3,6]] # only keep the parts relevant for 1D fit
            best_fit = None
            best_guess = np.inf
            # print('UPDATED GUESS: {}'.format(updated_guess))
            
            pUpper = [
                np.inf,
                50.0,
                np.max(r[0]),
                len(r[0]),
                40,
            ]
            pLower = [
                -np.inf,
                0.0,
                np.min(r[0]),
                0,
                -40,
            ]
            try:
                # print("Trying guess: {}".format(guess))
                res = least_squares(
                    gaussian_mask_sigma_1D,
                    updated_guess,
                    args=(np.array(r[0]), od_int),
                    bounds=(pLower, pUpper),
                    kwargs={
                        "mask_radius_x": EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                        "x0_mask": np.abs(updated_guess[2]),
                    },
                )
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
                resplot = resLSQ.x.copy()
                resplot[0] = 0 # set offset to zero for plotting, we add it later
                self.fittedImage = gaussian_mask_sigma_1D(
                    resplot,
                    np.array(r[0]),
                    0,
                    mask_radius_x=EXCLUSION_RADIUS * np.abs(updated_guess[3]),
                    x0_mask=np.abs(updated_guess[2]),
                ).reshape(od_int.shape) + resLSQ.x[0]
                # For plotting:
                fittedImage_no_mask = gaussian_mask_sigma_1D(
                    resplot,
                    np.array(r[0]),
                    0,
                    mask_radius_x=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                ).reshape(od_int.shape) + resLSQ.x[0]

            # Integration of number from computed column density

            # Compute the number by summing the pixels and multiplying by the pixel area
            # subtract the offset from the initial Gaussian fit, as well as the gradients
            # We make a background image where we just set the amplitude to zero, so that
            # we only have the offset and the gradients
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
            bg_res = [
                res2D[0],
                0,
                res2D[2],
                res2D[3],
                res2D[4],
                res2D[5],
                res2D[6],
                res2D[7],
            ]

            # fitted bg converted to density
            getsigma = SIGMA_0[self.odImage.species] * (
                2 if IMFIT_MODES[self.odImage.mode]["Image Path"] == "Vertical" else 1
            )

            fitted_bg = (
                gaussian_mask_sigma_no_rot(
                    bg_res,
                    r,
                    0,
                    mask_radius_x=1e-9,
                    mask_radius_y=1e-9,
                    x0_mask=np.abs(updated_guess[2]),
                    y0_mask=np.abs(updated_guess[4]),
                ).reshape(self.odImage.ODCorrected.shape)
                / getsigma
            )

            raw_number = (self.odImage.n - fitted_bg).sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            error = (
                np.sqrt((self.odImage.nerr**2).sum())
                * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            self.fitData = np.append(self.fitData, number)
            self.fitData = np.append(self.fitData, error)
            self.fitData = np.append(self.fitData, EXCLUSION_RADIUS)

            ### Get radial average
            print('FITDATA: {}'.format(self.fitData))

            # self.slices.radSlice = azimuthalAverage(self.odImage.ODCorrected, center)
            # self.slices.radSliceFit = azimuthalAverage(self.fittedImage, center)
            # self.slices.radSliceFitGauss = [0] * len(self.slices.radSlice)

            ### Calculate slices through fit
            y_for_x_slice = res2D[4] - r[1][0]
            # x_for_y_slice = res2D[2] - r[0][0]

            self.slices.points0 = od_int
            self.slices.ch0 = np.ones(len(r[0])) * y_for_x_slice + r[1][0]
            self.slices.fit0 = self.fittedImage
            self.slices.fit0a = fittedImage_no_mask

            # The other slice is just the cut along x as well
            self.slices.points1 = self.odImage.ODCorrected[
                round(y_for_x_slice), :
            ]
            self.slices.ch1 = None

            self.fitDataConfGauss = confidenceIntervals(resLSQ)
            self.fitDataGauss = resLSQ.x
            self.fittedImageGauss = gaussian1D(resLSQ.x, np.array(r[0]), 0)

            print("Done with fit function!")

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

        elif self.fitFunction == FIT_FUNCTIONS.index("Integrate"):
            # Integration of number from computed column density

            # Calculate average number density in border and subtract from rest of image
            border = int(max(min(self.odImage.n.shape) / 10, 5))
            border_mask = np.ones(self.odImage.n.shape)
            border_mask[border:-border, border:-border] = 0
            offset = np.sum(self.odImage.n * border_mask) / np.sum(border_mask)
            self.odImage.n -= offset

            interior = self.odImage.n[border:-border, border:-border]
            interior_err = self.odImage.nerr[border:-border, border:-border]

            # Compute the number by summing the pixels and multiplying by the pixel area

            raw_number = interior.sum()
            number = (
                raw_number * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            error = (
                np.sqrt((interior_err**2).sum())
                * (self.config["Pixel Size"] * self.odImage.data.bin) ** 2
            )

            # Compute the central position and size
            xrange = np.arange(interior.shape[1])
            yrange = np.arange(interior.shape[0])
            xc = np.sum(interior.sum(axis=0) * xrange) / raw_number
            yc = np.sum(interior.sum(axis=1) * yrange) / raw_number

            # Box size for plotting
            self.box = [
                [
                    self.odImage.xRange0.start + border,
                    self.odImage.xRange0.stop - border,
                ],
                [
                    self.odImage.xRange1.start + border,
                    self.odImage.xRange1.stop - border,
                ],
            ]

            # These calculations of moments aren't exact, since the contents of the square root can be negative since the column density can be negative in some pixels.
            # Therefore, all pixels less than 10 percent of the peak are dropped.
            sumx = interior.sum(axis=0)
            sumy = interior.sum(axis=1)
            minx = maximum(sumx) / 10
            miny = maximum(sumy) / 10
            sigx = np.sqrt(
                np.sum((sumx * (xrange - xc) ** 2)[sumx > minx])
                / np.sum(sumx[sumx > minx])
            )
            sigy = np.sqrt(
                np.sum((sumy * (yrange - yc) ** 2)[sumy > miny])
                / np.sum(sumy[sumy > miny])
            )

            xc += border + self.odImage.xRange0.start
            yc += border + self.odImage.xRange1.start

            xc = min(max(self.odImage.xRange0.start, xc), self.odImage.xRange0.stop - 1)
            yc = min(max(self.odImage.xRange1.start, yc), self.odImage.xRange1.stop - 1)

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

            self.slices.points0 = self.odImage.n[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = None

            self.slices.points1 = self.odImage.n[:, I0]
            self.slices.ch1 = [self.odImage.xRange0[I0]] * len(self.odImage.xRange1)
            self.slices.fit1 = None
        
        elif self.fitFunction == FIT_FUNCTIONS.index("Azimuthal Average Gauss"):
            # Integrate azimuthally. We first fit a 2D Gaussian, then we average over
            # an elliptical contour as determined from the Gaussian fit

            r = [None, None]
            r[0] = self.odImage.xRange0
            r[1] = self.odImage.xRange1

            print((r[1].start + r[1].stop) // 2)

            # INITIAL Gaussian fit with gradient
            ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
            p0 = [
                0,
                1,
                (r[0].start + r[0].stop) // 2,
                10,
                (r[1].start + r[1].stop) // 2,
                10,
                0,
                0,
                0
            ]

            # limit between 0 and 45 degrees
            pLower = [
                -np.inf,
                0,
                0,
                0,
                0,
                0,
                0,
                -np.inf,
                -np.inf,
            ]

            pUpper = [
                np.inf,
                np.inf,
                np.max(r[0]),
                len(r[0]),
                np.max(r[1]),
                len(r[1]),
                np.pi/4,
                np.inf,
                np.inf,
            ]

            resLSQ_G = least_squares(
                gaussianGradient,
                p0,
                args=(r, self.odImage.ODCorrected),
                bounds=(pLower, pUpper),
            )

            fit_res = gaussianGradient(resLSQ_G.x, r, 0).reshape(np.shape(self.odImage.ODCorrected))

            ### Get radial average
            x0f = resLSQ_G.x[2]
            y0f = resLSQ_G.x[4]
            sigx = resLSQ_G.x[3]
            sigy = resLSQ_G.x[5]
            theta = resLSQ_G.x[6]
                       
            I0 = self.odImage.xRange0.index(round(x0f))
            I1 = self.odImage.xRange1.index(round(y0f))

            X0, Y0 = np.meshgrid(r[0], r[1])
            XR = X0 * np.cos(theta) - Y0 * np.sin(theta)
            YR = X0 * np.sin(theta) + Y0 * np.cos(theta)

            x0R = x0f * np.cos(theta) - y0f * np.sin(theta)
            y0R = x0f * np.sin(theta) + y0f * np.cos(theta)
            R = np.sqrt((XR-x0R)**2/sigx**2 + (YR-y0R)**2/sigy**2) * np.sqrt(sigx*sigy)

            # Azimuthal averaging - the built in function does not allow elliptical averaging
            # calculate the mean - even if part of the annulus partially goes out of bounds
            # you're calculating the mean, so it's just averaged over less pixels
            f = lambda r : self.odImage.ODCorrected[(R >= r-.5) & (R < r+.5)].mean()
            ffit = lambda r : fit_res[(R >= r-.5) & (R < r+.5)].mean()
            
            r_az  = np.linspace(1, np.floor(np.max(R)) - 1)

            center = [I0, I1]

            ### Parameters: [offset, amplitude, wx, dODdx]
            p0 = [
                0,
                resLSQ_G.x[1],
                np.sqrt(sigx * sigy),
                0,
            ]

            pLower = [
                -np.inf,
                0,
                0,
                -np.inf,
            ]

            pUpper = [
                np.inf,
                np.inf,
                np.max(r[0]),
                np.inf,
            ]
            
            resLSQ = least_squares(
                azimGauss,
                p0,
                args=(r_az, np.vectorize(f)(r_az)),
                bounds=(pLower, pUpper),
            )
            
            offset = resLSQ.x[0]
            amplitude = resLSQ.x[1]
            sigx = resLSQ.x[2]
            dODdx = resLSQ.x[3]

            self.fitData = [offset, amplitude, sigx, dODdx]
            azimGaussfit = azimGauss(resLSQ.x, r_az, 0)

            self.slices.radSlice = np.vectorize(f)(r_az)
            self.slices.radSliceFit = azimGaussfit
            self.slices.radSliceFitGauss = np.vectorize(ffit)(r_az)

            ### Calculate slices through fit
            self.slices.points0 = self.odImage.ODCorrected[I1, :]
            self.slices.ch0 = [self.odImage.xRange1[I1]] * len(self.odImage.xRange0)
            self.slices.fit0 = fit_res[I1, :]

            self.slices.points1 = None
            self.slices.ch1 = None
            self.slices.fit1 = None

        else:
            print("Fit function undefined! Something went wrong!")
            return -1


class processFitResult:

    def __init__(self, fitObject, mode):

        self.fitObject = fitObject
        self.bin = self.fitObject.odImage.data.bin
        # self.atom = self.fitObject.odImage.atom
        self.config = IMFIT_MODES[mode]
        self.pixelSize = self.config["Pixel Size"]
        self.angle = self.config["Fit angle"]

        self.data = None
        self.data_dict = None

        self.getResults()

    def getResults(self):
        if self.fitObject.fitFunction == FIT_FUNCTIONS.index("Rotated Gaussian"):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": self.fitObject.fitData[6] * 180.0 / np.pi,
                "dODdx": 0,
                "dODdy": 0,
            }

            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
            ]
            self.data_dict = r

            rErr = {
                "offset": self.fitObject.fitDataConf[0],
                "peakOD": self.fitObject.fitDataConf[1],
                "x0": self.fitObject.fitDataConf[2],
                "y0": self.fitObject.fitDataConf[4],
                "wx": self.fitObject.fitDataConf[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitDataConf[5] * self.bin * self.pixelSize,
                "angle": self.fitObject.fitDataConf[6] * 180.0 / np.pi,
            }

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Gaussian"):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": 0,
                "dODdx": 0,
                "dODdy": 0,
            }

            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Gaussian Mask Sigma"):

            r = {
                "offset": self.fitObject.fitData[0],
                "dODdx": self.fitObject.fitData[7] / self.bin * self.pixelSize,
                "dODdy": self.fitObject.fitData[8] / self.bin * self.pixelSize,
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": self.fitObject.fitData[6],
                "N": self.fitObject.fitData[9],
                "Nerr": self.fitObject.fitData[10],
                "ExclR": self.fitObject.fitData[11],
            }

            print(r)
            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
                r["N"],
                r["Nerr"],
                r["ExclR"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index(
            "Gaussian Mask Sigma No Rot"
        ):

            r = {
                "offset": self.fitObject.fitData[0],
                "dODdx": self.fitObject.fitData[6] / self.bin * self.pixelSize,
                "dODdy": self.fitObject.fitData[7] / self.bin * self.pixelSize,
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": 0,
                "N": self.fitObject.fitData[8],
                "Nerr": self.fitObject.fitData[9],
                "ExclR": self.fitObject.fitData[10],
            }

            print(r)
            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
                r["N"],
                r["Nerr"],
                r["ExclR"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index(
            "Gaussian w/ Gradient"
        ) or self.fitObject.fitFunction == FIT_FUNCTIONS.index("Gaussian Fixed"):

            r = {
                "offset": self.fitObject.fitData[0],
                "dODdx": self.fitObject.fitData[6] / self.bin * self.pixelSize,
                "dODdy": self.fitObject.fitData[7] / self.bin * self.pixelSize,
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": 0,
            }

            print(r)
            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Twisted Gaussian"):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "angle": self.angle,
                "dODdx": 0,
                "dODdy": 0,
            }

            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Bigaussian"):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakODBEC": self.fitObject.fitData[1],
                "wxBEC": self.fitObject.fitData[2] * self.bin * self.pixelSize,
                "wyBEC": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "peakODThermal": self.fitObject.fitData[4],
                "wxThermal": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "wyThermal": self.fitObject.fitData[6] * self.bin * self.pixelSize,
                "x0": self.fitObject.fitData[7],
                "y0": self.fitObject.fitData[8],
            }

            self.data = [
                "fileName",
                r["peakODBEC"],
                r["wxBEC"],
                r["wyBEC"],
                r["peakODThermal"],
                r["wxThermal"],
                r["wyThermal"],
                r["x0"],
                r["y0"],
                r["offset"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac"):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "q": self.fitObject.fitData[6],
                "TTF": getTTF(self.fitObject)[0][0],
                "wxClassical": self.fitObject.fitDataGauss[3]
                * self.bin
                * self.pixelSize,
                "wyClassical": self.fitObject.fitDataGauss[5]
                * self.bin
                * self.pixelSize,
                "peakODClassical": self.fitObject.fitDataGauss[1],
            }

            self.data = [
                "fileName",
                r["peakODClassical"],
                r["wxClassical"],
                r["wyClassical"],
                r["peakOD"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["TTF"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index(
            "Fermi-Dirac fixed betamu"
        ):

            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "q": self.fitObject.fitData[6],
                "TTF": getTTF(self.fitObject)[0],
                "wxClassical": self.fitObject.fitDataGauss[3]
                * self.bin
                * self.pixelSize,
                "wyClassical": self.fitObject.fitDataGauss[5]
                * self.bin
                * self.pixelSize,
                "peakODClassical": self.fitObject.fitDataGauss[1],
                "TOF": self.fitObject.fitData[7],
                "fx": self.fitObject.fitData[8],
                "fy": self.fitObject.fitData[9],
                "fz": self.fitObject.fitData[10],
                "N": self.fitObject.fitData[11],
                "Nscaler": self.fitObject.fitData[12],
                "mass": self.fitObject.fitData[13],
            }

            self.data = [
                "fileName",
                r["peakODClassical"],
                r["wxClassical"],
                r["wyClassical"],
                r["peakOD"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["TTF"],
                r["TOF"],
                r["fx"],
                r["fy"],
                r["fz"],
                r["N"],
                r["Nscaler"],
                r["mass"],
            ]
            self.data_dict = r

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac 2D"):
            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "q": self.fitObject.fitData[6],
                "TTF": 0.5 * np.sqrt(np.pi / fermi_poly2(self.fitObject.fitData[6])[0]),
                "wxClassical": self.fitObject.fitDataGauss[3]
                * self.bin
                * self.pixelSize,
                "wyClassical": self.fitObject.fitDataGauss[5]
                * self.bin
                * self.pixelSize,
                "peakODClassical": self.fitObject.fitDataGauss[1],
            }

            self.data = [
                "fileName",
                r["peakODClassical"],
                r["wxClassical"],
                r["wyClassical"],
                r["peakOD"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["TTF"],
            ]
            self.data_dict = r
        
        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Fermi-Dirac 2D Int"):
            r = {
                "offset": self.fitObject.fitDataGauss[1],  # self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": 0,  # self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": 0,  # self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "q": self.fitObject.fitData[4],
                "TTF": np.sqrt(
                    -1 / (2 * mp.fp.polylog(2, -np.exp(self.fitObject.fitData[4])))
                ),
                "wxClassical": self.fitObject.fitDataGauss[3]
                * self.bin
                * self.pixelSize,
                "wyClassical": 0  # self.fitObject.fitDataGauss[5]
                * self.bin
                * self.pixelSize,
                "peakODClassical": self.fitObject.fitDataGauss[1],
                "gradient": self.fitObject.fitDataGauss[4],
            }

            self.data = [
                "fileName",
                r["peakODClassical"],
                r["wxClassical"],
                r["wyClassical"],
                r["peakOD"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["TTF"],
            ]
            self.data_dict = r
        
        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Twisted Fermi-Dirac 2D Int"):
            r = {
                "offset": self.fitObject.fitDataGauss[0],
                "peakOD": self.fitObject.fitDataGauss[1],
                "x0cl": self.fitObject.fitDataGauss[2],
                "wxcl": self.fitObject.fitDataGauss[3] * self.bin * self.pixelSize,
                "y0cl": self.fitObject.fitDataGauss[4],
                "wycl": self.fitObject.fitDataGauss[5] * self.bin * self.pixelSize,
                "offset_x": self.fitObject.fitData[0],
                "peakOD_x": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "q_x": self.fitObject.fitData[4],
                "TTF_x": np.sqrt(
                    -1 / (2 * mp.fp.polylog(2, -np.exp(self.fitObject.fitData[4])))
                ),
                "offset_y": self.fitObject.fitData[6],
                "peakOD_y": self.fitObject.fitData[7],
                "y0": self.fitObject.fitData[8],
                "wy": self.fitObject.fitData[9] * self.bin * self.pixelSize,
                "q_y": self.fitObject.fitData[10],
                "TTF_y": np.sqrt(
                    -1 / (2 * mp.fp.polylog(2, -np.exp(self.fitObject.fitData[10])))
                ),
                "number": self.fitObject.fitData[12],
            }

            self.data = [
                "fileName",
                r["peakOD"],
                r["wxcl"],
                r["wycl"],
                r["x0cl"],
                r["y0cl"],
                r["peakOD_x"],
                r["peakOD_y"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset_x"],
                r["offset_y"],
                r["TTF_x"],
                r["TTF_y"],
                r["number"],
            ]
            self.data_dict = r


        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Gauss (Mask) Int"):
            r = {
                "offset": self.fitObject.fitData[0],
                "dODdx": self.fitObject.fitData[4] / self.bin * self.pixelSize,
                "dODdy": 0,
                "peakOD": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": 0,
                "wx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "wy": 0,
                "angle": 0,
                "N": self.fitObject.fitData[5],
                "Nerr": self.fitObject.fitData[6],
                "ExclR": self.fitObject.fitData[7],
            }

            print(r)
            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["dODdy"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
                r["angle"],
                r["N"],
                r["Nerr"],
                r["ExclR"],
            ]
            self.data_dict = r


        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Thomas-Fermi"):
            r = {
                "offset": self.fitObject.fitData[0],
                "peakODTF": self.fitObject.fitData[1],
                "x0": self.fitObject.fitData[2],
                "y0": self.fitObject.fitData[4],
                "rx": self.fitObject.fitData[3] * self.bin * self.pixelSize,
                "ry": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "peakODGauss": self.fitObject.fitData[6],
                "sigxGauss": self.fitObject.fitData[7] * self.bin * self.pixelSize,
                "sigyGauss": self.fitObject.fitData[8] * self.bin * self.pixelSize,
            }

            self.data = [
                "fileName",
                r["peakODTF"],
                r["rx"],
                r["ry"],
                r["peakODGauss"],
                r["sigxGauss"],
                r["sigyGauss"],
                r["x0"],
                r["y0"],
                r["offset"],
            ]
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

        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Integrate"):
            r = {
                "offset": self.fitObject.fitData[0],
                "number": self.fitObject.fitData[1],
                "error": self.fitObject.fitData[2],
                "x0": self.fitObject.fitData[3],
                "y0": self.fitObject.fitData[4],
                "wx": self.fitObject.fitData[5] * self.bin * self.pixelSize,
                "wy": self.fitObject.fitData[6] * self.bin * self.pixelSize,
            }

            self.data = [
                "fileName",
                r["number"],
                r["error"],
                r["wx"],
                r["wy"],
                r["x0"],
                r["y0"],
                r["offset"],
            ]
            self.data_dict = r
        
        elif self.fitObject.fitFunction == FIT_FUNCTIONS.index("Azimuthal Average Gauss"):
            r = {
                "offset": self.fitObject.fitData[0],
                "peakOD": self.fitObject.fitData[1],
                "wx": self.fitObject.fitData[2] * self.bin * self.pixelSize,
                "dODdx": self.fitObject.fitData[3] / self.bin * self.pixelSize,
            }

            self.data = [
                "fileName",
                r["peakOD"],
                r["dODdx"],
                r["wx"],
                r["offset"],
            ]
            self.data_dict = r

        else:
            print("Fit function undefined! Something went wrong!")
            return -1

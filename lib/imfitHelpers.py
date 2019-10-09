from imfitDefaults import *
import polylog
import numpy as np

def checkAtom(atom):

    if isinstance(atom, basestring):
        atom = atom.upper()
        if atom not in ATOM_NAMES:
            print('Error: Atom must be one of: ' + ', '.join(ATOM_NAMES))
            return -1
        else:
            return ATOM_NAMES.index(atom)

    elif isinstance(atom, int):
        if atom not in range(NATOMS):
            print('Error: Atom index must be in the range 0-{0}'.format(NATOMS-1))
            return -1
        else:
            return atom

    else:
        print('Atom must be a string or integer corresponding to an atom in the experiment.')
        return -1

def checkFrame(frame):

        if isinstance(frame, basestring):
            frame = frame.lower()
            if frame not in FRAMESEQUENCE:
                print('Error: Frame must be one of: ' + ', '.join(FRAMESEQUENCE))
                return -1
            else:
                return FRAMESEQUENCE.index(frame)
        elif isinstance(frame, int):
            if frame not in range(self.nFrames):
                print('Error: Frame index must be in the range 0-{0}'.format(self.nFrames))
                return -1
            else:
                return frame
        else:
            print('Frame must be a string or integer corresponding to a frame in the experiment.')

def checkFitFunction(fitFunction):

    if isinstance(fitFunction, basestring):
        if fitFunction in FIT_FUNCTIONS:
            return FIT_FUNCTIONS.index(fitFunction)
        else:
            print("Error: Fit function must be one of " ','.join(FIT_FUNCTIONS))
            return -1
    elif isinstance(fitFunction, int):
        if fitFunction < len(FIT_FUNCTIONS):
            return fitFunction
        else:
            print("Error: only {} fitting functions are defined!".format(len(FIT_FUNCTIONS)))
            return -1
    else:
        print("Error: Fit function must be a string or integer corresponding to a real fitting function.")
        return -1

def cropArray(array, x, y):
    ### Crop array according to values in a list
    t = array[x,:]
    return t[:,y]

def confidenceIntervals(res_lsq):
    hess = np.matmul(res_lsq.jac.transpose(),res_lsq.jac)
    try:
        cov = np.linalg.inv(hess) * (res_lsq.fun**2).mean()
        return np.sqrt(np.abs(np.diag(cov)))

    except np.linalg.LinAlgError as err:
        if "Singular matrix" in str(err):
            print('Error: Least squares method returned a singular hessian due to a poor fit. Fit uncertainty not estimated. Zeros returned.')
            return np.array([0]*hess.shape[0])
        else:
            return -1
    


def getTTF(fitObject):
    if fitObject.fitFunction == FIT_FUNCTIONS.index('Fermi-Dirac'):
        TTF = (6.0 * polylog.fermi_poly3(fitObject.fitData[6]))**(-1.0/3.0)
        TTFErr = 0.5*( (6.0 * polylog.fermi_poly3(fitObject.fitData[6]-fitObject.fitDataConf[6]))**(-1.0/3.0) - (6.0 * polylog.fermi_poly3(fitObject.fitData[6]+fitObject.fitDataConf[6]))**(-1.0/3.0))
        return TTF,TTFErr
    else:
        print('T/TF only available for Fermi-Dirac fit.')
        return -1, -1

def azimuthalAverage(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x-center[0])**2.0 + (y-center[1])**2.0)
    r = r.astype(np.int)
    
    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())

    return tbin/nr


def getLastFile(path):

    import os
    d = os.listdir(path)

    d1 = [i.split('.')[0] for i in d]
    d2 = [i.split('_')[-1] for i in d1]

    if not d2:
        return 0
    else:
        f = max(map(int, d2))
        return f

def getImagesFromRange(stringIn):

    l = []
    
    l1 = stringIn.split(',')
    for k in l1:
        if '-' in k:
            q = k.split('-')
            for m in range(int(q[0]),int(q[1])+1):
                l.append(int(m))
        elif ':' in k:
            q = k.split(':')
            for m in range(int(q[0]),int(q[1])+1):
                l.append(int(m))
        else:
            l.append(int(k))

    l.sort()
    
    return l

def upload2Origin(atom, fitFunction, data):
    import win32com.client

    progID = 'Origin.ApplicationSI'
    orgApp = win32com.client.Dispatch(progID)


    atom = checkAtom(atom)
    fitFunction = checkFitFunction(fitFunction)



    if atom != -1 and fitFunction != -1:

        if atom == 2:
            if FIT_FUNCTIONS[fitFunction] == 'Gaussian':
                template = 'KRbGauss1'
            elif FIT_FUNCTIONS[fitFunction] == 'Twisted Gaussian':
                template = 'KRbGauss1'
            else:
                template = WORKSHEET_NAMES[fitFunction]    
        else:
            template = WORKSHEET_NAMES[fitFunction]
        worksheetName = ATOM_NAMES[atom] + WORKSHEET_NAMES[fitFunction]

        longname = ATOM_NAMES[atom] + " " + FIT_FUNCTIONS[fitFunction]

        if orgApp.FindWorksheet(worksheetName) is None:
            orgApp.CreatePage(2, worksheetName, template)
        orgApp.Execute("{}!page.longname$ = {}".format(worksheetName, longname))
        # orgApp.Execute("{}!page.active$ = {}".format(worksheetName, "Sheet1")) 


        n = 0
        for k in data:
            uploadSuccess = orgApp.PutWorksheet("[{}]Sheet1".format(worksheetName), k, -1, n)
            if uploadSuccess:
                n += 1
            else:
                print("Failed to upload to Origin. Is Sheet1 in the proper workbook available?")
                return -1

    else: 
        print("Failed to upload to Origin.")
        return -1

def checkGuess(x0, xUpper, xLower):
    n = len(x0)
    for k in range(n):

        if np.isnan(x0[k]):
            x0[k] = xUpper[k]

        if x0[k] > xUpper[k]:
            x0[k] = xUpper[k]
            print("Initial guess clamped to upper bound.")

        elif x0[k] < xLower[k]:
            x0[k] = xLower[k]
            print("Initial guess clamped to lower bound.")
        
        else:
            pass
    
    return x0

if __name__ == "__main__":
    uploadToOrigin('K', 'Gaussian', [1,1])

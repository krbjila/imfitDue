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
        return -2

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

def cropArray(array, x, y):
    ### Crop array according to values in a list
    t = array[x,:]
    return t[:,y]

def confidenceIntervals(res_lsq):
    cov = np.linalg.inv(np.matmul(res_lsq.jac.transpose(),res_lsq.jac)) * (res_lsq.fun**2).mean()
    return np.sqrt(np.abs(np.diag(cov)))


def getTTF(fitObject):
    if FIT_FUNCTIONS[fitObject.fitFunction] == 'fermidirac':
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

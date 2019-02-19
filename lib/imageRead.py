from imfitDefaults import *
from imfitHelpers import *
from spe import SpeFile
import numpy as np
import os 

### Define separate classes for reading in images from either the Princeton Instruments
### camera or the Andor iXon 888. Regardless of the camera, the classes have the same
### member variables, and the data is structured in the same way.

class readPIImage():

    def __init__(self,path):
        self.path = path

        self.nFrames = None
        self.hImgSize = None
        self.vImgSize = None
        self.bin = None

        self.fileLocation = None
        self.fileName = None 
        self.fileNumber = None
        self.fileTimeStamp = None

        self.camera = CAMERA_NAME_PI
        self.img = None

        self.updateAll()

    def setPath(self,path):
        self.path = path 
        return self.updateAll()

    def updateAll(self):
    
        if not os.path.isfile(self.path):
            print('No such file or directory: ' + self.path)
            return 0 

        self.fileLocation = FILESEP.join(self.path.split(FILESEP)[0:-1]) + FILESEP,
        self.fileName = self.path.split(FILESEP)[-1],
        self.fileNumber =  int(self.path.split('_')[-1].split('.')[0]),
        self.fileTimeStamp = os.path.getctime(self.path)

        x = SpeFile(self.path)
        
        self.nFrames = int(x.data.shape[0])
        self.hImgSize = int(x.data.shape[2]/NATOMS) #Account for the number of shifts
        self.vImgSize = int(x.data.shape[1])

        if self.hImgSize*NATOMS == SENSOR_WIDTH_PI:
            self.bin = False
        elif self.hImgSize*NATOMS == SENSOR_WIDTH_PI/2:
            self.bin = True
        else:
            print('Oh No! Something weird happened with the binning of pixels!')


        self.img = np.zeros((int(NATOMS*self.vImgSize*self.nFrames),int(self.hImgSize)))

        for n in range(NATOMS):
            for k in range(self.nFrames):
                self.img[(k + n*self.nFrames)*self.vImgSize:(k + n*self.nFrames+1)*self.vImgSize,:] = x.data[k,:,n*self.hImgSize:(n+1)*self.hImgSize]

        return 1

    def getFrame(self, atom, frame):
        atom = checkAtom(atom)
        frame = checkFrame(frame)

        return self.img[(frame + atom*self.nFrames)*self.vImgSize:(frame + atom*self.nFrames+1)*self.vImgSize,:]

class readIXONImage():

    def __init__(self, path):
        self.path = path

        self.nFrames = None
        self.hImgSize = None
        self.vImgSize = None

        self.fileLocation = None
        self.fileName = None 
        self.fileNumber = None
        self.fileTimeStamp = None

        self.camera = CAMERA_NAME_IXON
        self.img = None

        self.updateAll()

    def setPath(self,path):
        self.path = path 
        return self.updateAll()

    def updateAll(self):
    
        if not os.path.isfile(self.path):
            print('No such file or directory: ' + self.path)
            return 0 

        self.fileLocation = FILESEP.join(self.path.split(FILESEP)[0:-1]) + FILESEP,
        self.fileName = self.path.split(FILESEP)[-1],
        self.fileNumber =  int(self.path.split('_')[-1].split('.')[0]),
        self.fileTimeStamp = os.path.getctime(self.path)

        self.img = np.genfromtxt(self.path, delimiter=',')
        
        self.nFrames = 3
        self.hImgSize = int(self.img.shape[1])
        self.vImgSize = int(self.img.shape[0]/(NATOMS*self.nFrames))

        if self.hImgSize == SENSOR_WIDTH_IXON:
            self.bin = False
        elif self.hImgSize == SENSOR_WIDTH_IXON/2:
            self.bin = True
        else:
            print('Oh No! Something weird happened with the binning of pixels!')

        return 1

    def getFrame(self, atom, frame):
        atom = checkAtom(atom)
        frame = checkFrame(frame)

        return self.img[(frame + atom*self.nFrames)*self.vImgSize:(frame + atom*self.nFrames+1)*self.vImgSize,:]


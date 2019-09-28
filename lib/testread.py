from imfitDefaults import *
from imfitHelpers import *
from spe import SpeFile
import numpy as np
import os 


class readXimeaImage():

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

        self.camera = CAMERA_NAME_XIMEA
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

        import json 

        f = open(self.path,'r')
        dataTemp = json.load(f)
        
        

        self.img = dataTemp['K']['Shadow']
        self.img = np.vstack((self.img, dataTemp['K']['Bright']))
        self.img = np.vstack((self.img, dataTemp['K']['Dark']))
        self.img = np.vstack((self.img, dataTemp['Rb']['Shadow']))
        self.img = np.vstack((self.img, dataTemp['Rb']['Bright']))
        self.img = np.vstack((self.img, dataTemp['Rb']['Dark']))

        # self.img = np.genfromtxt(self.path, delimiter=',')
        
        self.nFrames = 3
        self.hImgSize = int(self.img.shape[1])
        self.vImgSize = int(self.img.shape[0]/(NATOMS*self.nFrames))


        self.bin = False        
        # if self.hImgSize == SENSOR_WIDTH_IXON:
        #     self.bin = False
        # elif self.hImgSize == SENSOR_WIDTH_IXON/2:
        #     self.bin = True
        # else:
        #     print('Oh No! Something weird happened with the binning of pixels!')

        return 1

    def getFrame(self, atom, frame):
        atom = checkAtom(atom)
        frame = checkFrame(frame)

        return self.img[(frame + atom*self.nFrames)*self.vImgSize:(frame + atom*self.nFrames+1)*self.vImgSize,:]



path = "K:/data/2019/09/20190926/ximea/XI20190926_004.dat"

readXimeaImage(path)
from imfitDefaults import *
from imfitHelpers import *
from spe import SpeFile
import numpy as np
import os 

### Define separate classes for reading in images from either the Princeton Instruments
### camera or the Andor iXon 888. Regardless of the camera, the classes have the same
### member variables, and the data is structured in the same way.

def readImage(mode, path):
    extension = IMFIT_MODES[mode]['Extension Filter']

    if extension == '*.csv':
        reader = CsvReader(mode)
    elif extension == '*.npz':
        reader = NpzReader(mode)
    reader.setPath(path)
    return reader

class Reader(object):
    def __init__(self, mode):
        self.mode = mode

        self.config = IMFIT_MODES[self.mode]
        self.path = self.config['Default Path']
        self.metadata = {}
        self.frames = {}

        self.n_frames = self.config['Number of Frames']
        self.frame_order = self.config['Frame Order']

        # Need to know the size of the array if we are to auto detect binning
        self.auto_detect_binning = self.config['Auto Detect Binning']
        if self.auto_detect_binning:
            self.array_width = self.config['Array Width']
        
        # self.nFrames = None
        # self.hImgSize = None
        # self.vImgSize = None
        # self.bin = None

    def setPath(self,path):
        self.path = path 
        return self.updateAll()

    def updateAll(self):
        if not os.path.isfile(self.path):
            print('No such file or directory: ' + self.path)
            return 0

        self.fileLocation = FILESEP.join(self.path.split(FILESEP)[0:-1]) + FILESEP
        self.fileName = self.path.split(FILESEP)[-1]

        # Convention: "filebase_shot.extension"
        self.fileNumber =  int(self.path.split('_')[-1].split('.')[0])
        self.fileTimeStamp = os.path.getctime(self.path)
        
        (frame_list, metadata) = self.getData()

        self.frames  = {}
        for species, frame_dict in self.config['Frame Order'].items():
            temp = {}
            for frame_name, frame_index in frame_dict.items():
                temp.update({frame_name: np.array(frame_list[frame_index])})
            self.frames.update({species: temp})
        
        self.metadata = metadata
        self.bin = metadata['bins']
        self.hImgSize = metadata['hImgSize']
        self.vImgSize = metadata['vImgSize']

    def getFrame(self, species, frame):
        return self.frames[species][frame]

class CsvReader(Reader):
    def __init__(self, mode, delimiter=','):
        super(CsvReader, self).__init__(mode)
        self.delimiter = delimiter

    def getData(self):
        data = np.genfromtxt(self.path, delimiter=self.delimiter)

        # We don't really have metadata from csv's, except for binning info
        metadata = {}
        if self.auto_detect_binning:
            bins = self.array_width / int(data.shape[1])
            metadata['bins'] = bins
        else:
            # Assume no binning for now
            # TODO: check this intelligently
            metadata['bins'] = 1
        frame_size = data.shape[0] / self.n_frames # concatenated along vertical dimension
        frames = [data[i*frame_size:(i+1)*frame_size,:] for i in range(self.n_frames)]    

        metadata['hImgSize'] = data.shape[1]
        metadata['vImgSize'] = frame_size
        return (frames, metadata)

class NpzReader(Reader):
    def __init__(self, mode):
        super(NpzReader, self).__init__(mode)

    def getData(self):
        data = np.load(self.path, allow_pickle=True)
        return (data['data'], data['meta'])
    
class DatReader(Reader):
    def __init__(self, mode):
        super(DatReader, self).__init__(mode)

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

class readIXONGSMImage():

    def __init__(self, path):
        self.path = path

        self.nFrames = None
        self.hImgSize = None
        self.vImgSize = None

        self.fileLocation = None
        self.fileName = None 
        self.fileNumber = None
        self.fileTimeStamp = None

        self.camera = CAMERA_NAME_IXON_GSM
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

        
        if self.hImgSize == SENSOR_WIDTH_IXON_GSM:
            self.bin = False
        elif self.hImgSize == SENSOR_WIDTH_IXON_GSM/2:
            self.bin = True
        else:
            print('Oh No! Something weird happened with the binning of pixels!')

        return 1

    def getFrame(self, atom, frame):
        atom = checkAtom(atom)
        frame = checkFrame(frame)

        return self.img[(frame + atom*self.nFrames)*self.vImgSize:(frame + atom*self.nFrames+1)*self.vImgSize,:]

class readIXONVImage():

    def __init__(self, path):
        self.path = path

        self.nFrames = None
        self.hImgSize = None
        self.vImgSize = None

        self.fileLocation = None
        self.fileName = None 
        self.fileNumber = None
        self.fileTimeStamp = None

        self.camera = CAMERA_NAME_IXONV
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

        
        if self.hImgSize == SENSOR_WIDTH_IXONV:
            self.bin = False
        elif self.hImgSize == SENSOR_WIDTH_IXONV/2:
            self.bin = True
        else:
            print('Oh No! Something weird happened with the binning of pixels!')

        return 1

    def getFrame(self, atom, frame):
        atom = checkAtom(atom)
        frame = checkFrame(frame)

        return self.img[(frame + atom*self.nFrames)*self.vImgSize:(frame + atom*self.nFrames+1)*self.vImgSize,:]

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
        f.close()       

        self.img = np.array(dataTemp['K']['Shadow'], dtype=np.float)
        self.img = np.vstack((self.img, np.array(dataTemp['K']['Bright'], dtype=np.float)))
        self.img = np.vstack((self.img, np.array(dataTemp['K']['Dark'], dtype=np.float)))
        self.img = np.vstack((self.img, np.array(dataTemp['Rb']['Shadow'], dtype=np.float)))
        self.img = np.vstack((self.img, np.array(dataTemp['Rb']['Bright'], dtype=np.float)))
        self.img = np.vstack((self.img, np.array(dataTemp['Rb']['Dark'], dtype=np.float)))

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
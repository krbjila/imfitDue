import warnings
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
    else:
        print("Could not load {}: extension {} invalid.".format(path, extension))
        return None
    try:
        reader.setPath(path)
    except Exception as e:
        print("Could not load {}: {}".format(path, e))
        return None
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
        self.updateAll()

    def updateAll(self):
        if not os.path.isfile(self.path):
            raise(IOError('No such file or directory: ' + self.path))

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
        frame_size = int(data.shape[0] / self.n_frames) # concatenated along vertical dimension
        frames = [data[i*frame_size:(i+1)*frame_size,:] for i in range(self.n_frames)]    

        metadata['hImgSize'] = data.shape[1]
        metadata['vImgSize'] = frame_size
        return (frames, metadata)

class NpzReader(Reader):
    def __init__(self, mode):
        super(NpzReader, self).__init__(mode)

    def getData(self):
        f = np.load(self.path, allow_pickle=True)
        data = f['data']
        metadata = f['meta'].item()
        # Assume symmetric binning
        metadata['bins'] = metadata['binning'][0]
        metadata['hImgSize'] = data.shape[2]
        metadata['vImgSize'] = data.shape[1]
        # TODO: restructure metadata??? see Reader.updateAll()
        return (data, metadata)
    
    # TODO: Implement this!
class DatReader(Reader):
    def __init__(self, mode):
        super(DatReader, self).__init__(mode)

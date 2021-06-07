import datetime

# from matplotlib.pyplot import autoscale
now = datetime.datetime.now()

DEFAULT_PATH = now.strftime('K:/data/%Y/%m/%Y%m%d/') # The dataserver is mounted to K: on the Windows computers.
FILESEP = '/'

AUTOSCALE_MIN = 2 # percentile
AUTOSCALE_HEADROOM = 1.1 # factor above max

# FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian', 'Bigaussian', 'Fermi-Dirac', 'Vertical BandMap']
FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian', 'Bigaussian', 'Fermi-Dirac']
KRB_FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian']
NONTWISTED_FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian']
WORKSHEET_NAMES = [ 'GaussGrad', 'Gauss1', 'Gauss1', 'Gauss1', 'Gauss2', 'FermiDirac', 'BandMapV']

DEFAULT_MODE = 'Axial iXon'

from collections import OrderedDict
IMFIT_MODES = OrderedDict([
    ('Axial iXon', {
        ### Mode 1 - Axial iXon
        'Default Path': DEFAULT_PATH + 'Andor/',
        'Default Suffix': 'ixon_{}.npz',
        'Pixel Size': 2.58,
        'Species': ['K', 'Rb'],
        'Image Path': 'Axial',
        'Default Region': [[150, 220, 250, 250], 
                  [150, 350, 300, 300]],
        'Extension Filter': '*.npz',
        'Fit Functions': FIT_FUNCTIONS,
        'Enforce same fit for both': False,
        'Auto Detect Binning': True,
        'Array Width': 512,
        'Number of Frames': 6,
        'Frame Order': {
            'K': {
                'Shadow': 0,
                'Light': 1,
                'Dark': 2,
            },
            'Rb': {
                'Shadow': 3,
                'Light': 4,
                'Dark': 5
            }
        },
        'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        'CSat': {'K': 19e3 / 4.0, 'Rb': 19e3 / 4.0}, # Unbinned effective C_sat
    }),
    ('Axial iXon Molecules ToF', {
        ### Mode 1.1 - Axial iXon
        'Default Path': DEFAULT_PATH + 'KRbFK/',
        'Default Suffix': 'krbfk_{}.npz',
        'Pixel Size': 2.58,
        'Species': ['|0,0>', '|1,0>'],
        'Image Path': 'Axial',
        'Default Region': [[135, 150, 100, 100], 
                  [135, 150, 100, 100]],
        'Extension Filter': '*.npz',
        'Fit Functions': NONTWISTED_FIT_FUNCTIONS,
        'Enforce same fit for both': True,
        'Auto Detect Binning': True,
        'Array Width': 512,
        'Number of Frames': 4,
        'Frame Order': {
            '|0,0>': {
                'Shadow': 0,
                'Light': 2,
                'Dark': 1,
            },
            '|1,0>': {
                'Shadow': 0,
                'Light': 2,
                'Dark': 1
            }
        },
        'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        'CSat': {'|0,0>': 19e3 / 4.0, '|1,0>': 19e3 / 4.0}, # Unbinned effective C_sat
        'Allow fit both states': True
    }),
    ('Axial iXon Molecules In Situ', {
        ### Mode 3 - Axial iXon Molecules
        'Default Path': DEFAULT_PATH + 'MoleculeInSituFK/',
        'Default Suffix': 'ixon_{}.npz',
        'Pixel Size': 2.58,
        'Species': ['|0,0>', '|1,0>'],
        'Image Path': 'Axial',
        'Default Region': [[268, 253, 50, 25], 
                  [268, 253, 50, 25]],
        'Extension Filter': '*.npz',
        'Fit Functions': KRB_FIT_FUNCTIONS,
        'Enforce same fit for both': True,
                'Auto Detect Binning': True,
        'Array Width': 512,
        'Number of Frames': 6,
        'Frame Order': {
            '|0,0>': {
                'Shadow': 0,
                'Light': 1,
                'Dark': 2,
            },
            '|1,0>': {
                'Shadow': 3,
                'Light': 4,
                'Dark': 5
            }
        },
        'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        'CSat': {'|0,0>': 19e3 / 4.0, '|1,0>': 19e3 / 4.0}, # Unbinned effective C_sat
    }),
    ('Pixelfly Test', {
        ### PLACEHOLDER
        ### TODO: put in finalized values
        'Default Path': DEFAULT_PATH + 'Pixelfly/',
        'Default Suffix': 'pixelfly_{}.npz',
        'Pixel Size': 2.58,
        'Species': ['K', 'Rb'],
        'Image Path': 'Axial',
        'Default Region': [[150, 220, 250, 250], 
                  [150, 350, 300, 300]],
        'Extension Filter': '*.npz',
        'Fit Functions': FIT_FUNCTIONS,
        'Enforce same fit for both': False,
        'Auto Detect Binning': False,
        'Number of Frames': 6,
        'Frame Order': {
            'K': {
                'Shadow': 0,
                'Light': 2,
                'Dark': 4,
            },
            'Rb': {
                'Shadow': 1,
                'Light': 3,
                'Dark': 5
            }
        },
        'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        'CSat': {'K': 19e3 / 4.0, 'Rb': 19e3 / 4.0}, # Unbinned effective C_sat
    }),
    ('Vertical iXon', {
        ### Mode 2 - Vertical iXon
        'Default Path': DEFAULT_PATH + 'Andor_Vertical/',
        'Default Suffix': 'twospecies_{}.npz',
        'Pixel Size': 0.956,
        'Species': ['K', 'Rb'],
        'Image Path': 'Vertical',
        'Default Region': [[125, 280, 250, 250], 
                  [125, 280, 250, 250]],
        'Extension Filter': '*.npz',
        'Fit Functions': FIT_FUNCTIONS,
        'Enforce same fit for both': False,
        'Auto Detect Binning': True,
        'Array Width': 512,
        'Number of Frames': 6,
        'Frame Order': {
            'K': {
                'Shadow': 0,
                'Light': 1,
                'Dark': 2,
            },
            'Rb': {
                'Shadow': 3,
                'Light': 4,
                'Dark': 5
            }
        },
        'Fit angle': 32.0, # Deg, Twisted Gaussian fit
        'CSat': {'K': 1400, 'Rb': 1400}, # Unbinned effective C_sat
    }),
    # ('Ximea', {
    #     ### Mode 0 - Ximea
    #     'Default Path': DEFAULT_PATH + 'ximea/',
    #     'Default Suffix': 'xi_{}.csv',
    #     'Pixel Size': 3.46,
    #     'Species': ['K', 'Rb'],
    #     'Image Path': 'Side',
    #     'Default Region': [[185, 260, 300, 300], 
    #               [185, 260, 300, 300]],
    #     'Extension Filter': '*.dat',
    #     'Fit Functions': FIT_FUNCTIONS,
    #     'Enforce same fit for both': False
    # }),
])
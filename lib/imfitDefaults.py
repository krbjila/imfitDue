
#######################################################################
########################    Imaging Parameters   ######################
#######################################################################

IMAGING_PATHS = ['Axial', 'Side', 'Vertical']
PIXEL_SIZES = [2.58, 3.46, 0.956]
CSATEFF = [19000, 1.0E6, 1400.0]

FRAMESEQUENCE = ['shadow', 'light', 'dark']
ATOM_NAMES = ['K','RB','KRB', 'KRBINSITU']
NATOMS = 2
WAVELENGTHS = [767.0E-9, 780.0E-9]


VERT_TRAP_ANGLE = 30.0 # Degrees

#######################################################################
########################  File System Parameters  #####################
#######################################################################

import datetime
now = datetime.datetime.now()

# with open('./lib/ip.txt') as f:
# 	ip_str = f.read(100)
# DEFAULT_PATH = now.strftime('//'+ip_str+'/krbdata/data/%Y/%m/%Y%m%d/') # PolarKRB's IP address
DEFAULT_PATH = now.strftime('K:/data/%Y/%m/%Y%m%d/') # The dataserver is mounted to K: on the Windows computers.


DEFAULT_PATH_PI = DEFAULT_PATH + 'ximea/'
DEFAULT_PATH_IXON = DEFAULT_PATH + 'Andor/'
DEFAULT_PATH_IXON_GSM = DEFAULT_PATH + 'MoleculeInSituFK/'
# DEFAULT_PATH_IXON_GSM = DEFAULT_PATH + 'Andor/'
DEFAULT_PATH_IXONV = DEFAULT_PATH + 'Andor_Vertical/'
FILESEP = '/'

#######################################################################
########################  Princeton Instruments  ######################
#######################################################################

CAMERA_NAME_PI = 'Princeton Instruments'




SENSOR_WIDTH_PI = 808

TPROBE_PI = [40.0, 40.0]
ODSAT_PI = [6.0, 6.0]
ISAT_FLUX_PI = [290.0, 260.0] #Saturation intensity flux counts/px/microsecond (unbinned)

#######################################################################
########################      Andor iXon 888     ######################
#######################################################################

CAMERA_NAME_IXON = 'Andor iXon 888'

SENSOR_WIDTH_IXON = 512 # Was 1024, but rotated camera orientation so the skinny side of the kinetics is horizontal

TPROBE_IXON = [40.0, 40.0]
ODSAT_IXON = [4.0, 4.0]
ISAT_FLUX_IXON = [275.0, 275.0]

#######################################################################
########################      Andor iXon 888 Molecules     ############
#######################################################################

CAMERA_NAME_IXON_GSM = 'Andor iXon 888 (Molecules)'

SENSOR_WIDTH_IXON_GSM = 512 # Was 1024, but rotated camera orientation so the skinny side of the kinetics is horizontal

TPROBE_IXON_GSM = [40.0, 40.0]
ODSAT_IXON_GSM = [4.0, 4.0]
ISAT_FLUX_IXON_GSM = [275.0, 275.0]


#######################################################################
###################      Andor iXon 888 Vertical   ####################
#######################################################################

CAMERA_NAME_IXONV = 'Andor iXon 888 (Vertical)'

SENSOR_WIDTH_IXONV = 512 # Was 1024, but rotated camera orientation so the skinny side of the kinetics is horizontal

TPROBE_IXONV = [40.0, 40.0]
ODSAT_IXONV = [4.0, 4.0]
ISAT_FLUX_IXONV = [275.0, 275.0]

#######################################################################
########################        Ximea xiQ        ######################
#######################################################################

CAMERA_NAME_XIMEA = 'Ximea xiQ'
ODSAT_XIMEA = [6.0, 6.0]


########################################################################
######################     Fitting Parameters     ######################
########################################################################

DEFAULT_REGION_XIMEA = [[185, 260, 300, 300], 
                  [185, 260, 300, 300]]
DEFAULT_REGION_IXON = [[150, 220, 250, 250], 
                  [150, 350, 300, 300]]
DEFAULT_REGION_IXON_GSM = [[280, 240, 60, 40], 
                  [280, 240, 60, 40]]
DEFAULT_REGION_IXONV = [[125, 280, 250, 250], 
                  [125, 280, 250, 250]]
DEFAULT_REGION = {
	'XIMEA' : DEFAULT_REGION_XIMEA,
	'IXON' : DEFAULT_REGION_IXON,
	'IXON_GSM' : DEFAULT_REGION_IXON_GSM,
	'IXONV' : DEFAULT_REGION_IXONV
}

FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian', 'Bigaussian', 'Fermi-Dirac', 'Vertical BandMap']
KRB_FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian']
WORKSHEET_NAMES = [ 'GaussGrad', 'Gauss1', 'Gauss1', 'Gauss1', 'Gauss2', 'FermiDirac', 'BandMapV']

DEFAULT_MODE = 'Axial iXon'

from collections import OrderedDict
IMFIT_MODES = OrderedDict([
    ('Axial iXon', {
        ### Mode 1 - Axial iXon
        'Default Path': DEFAULT_PATH + 'Andor/',
        'Default Suffix': 'ixon_{}.csv',
        'Pixel Size': 2.58,
        'Species': ['K', 'Rb'],
        'Image Path': 'Axial',
        'Default Region': [[150, 220, 250, 250], 
                  [150, 350, 300, 300]],
        'Extension Filter': '*.csv',
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
        'CSat': {'K': 19e3 / 4.0, 'Rb': 19e3 / 4.0}, # Unbinned effective C_sat
    }),
    ('Axial iXon Molecules', {
        ### Mode 3 - Axial iXon Molecules
        'Default Path': DEFAULT_PATH + 'MoleculeInSituFK/',
        'Default Suffix': 'ixon_{}.csv',
        'Pixel Size': 2.58,
        'Species': ['|0,0>', '|1,0>'],
        'Image Path': 'Axial',
        'Default Region': [[280, 240, 60, 40], 
                  [280, 240, 60, 40]],
        'Extension Filter': '*.csv',
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
        'CSat': {'|0,0>': 19e3 / 4.0, '|1,0>': 19e3 / 4.0}, # Unbinned effective C_sat
    }),
    ('Vertical iXon', {
        ### Mode 2 - Vertical iXon
        'Default Path': DEFAULT_PATH + 'Andor_Vertical/',
        'Default Suffix': 'twospecies_{}.csv',
        'Pixel Size': 0.956,
        'Species': ['K', 'Rb'],
        'Image Path': 'Vertical',
        'Default Region': [[125, 280, 250, 250], 
                  [125, 280, 250, 250]],
        'Extension Filter': '*.csv',
        'Fit Functions': FIT_FUNCTIONS,
        'Enforce same fit for both': False
    }),
    ('Ximea', {
        ### Mode 0 - Ximea
        'Default Path': DEFAULT_PATH + 'ximea/',
        'Default Suffix': 'xi_{}.csv',
        'Pixel Size': 3.46,
        'Species': ['K', 'Rb'],
        'Image Path': 'Side',
        'Default Region': [[185, 260, 300, 300], 
                  [185, 260, 300, 300]],
        'Extension Filter': '*.dat',
        'Fit Functions': FIT_FUNCTIONS,
        'Enforce same fit for both': False
    }),
])
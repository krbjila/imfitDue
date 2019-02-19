
FILESEP = '/'
DEFAULT_PATH = '/home/fermikrb/luigi/imfit/testfiles/pi_192.spe'

#######################################################################
########################    Imaging Parameters   ######################
#######################################################################

FRAMESEQUENCE = ['shadow', 'light', 'dark']
ATOM_NAMES = ['K','RB']
NATOMS = len(ATOM_NAMES)
WAVELENGTHS = [767.0E-9, 780.0E-9]


IMAGING_PATHS = ['axial', 'axial low mag.', 'vertical']
PIXEL_SIZES = [2.3, 7.5, 1.0]

#######################################################################
########################  Princeton Instruments  ######################
#######################################################################

CAMERA_NAME_PI = 'Princeton Instruments'


DEFAULT_PATH_PI = '/home/fermikrb/luigi/imfit/testfiles/'

SENSOR_WIDTH_PI = 808

TPROBE_PI = [40.0, 40.0]
ODSAT_PI = [4.0, 4.0]
ISAT_FLUX_PI = [290.0, 260.0] #Saturation intensity flux counts/px/microsecond (unbinned)

#######################################################################
########################      Andor iXon 888     ######################
#######################################################################

CAMERA_NAME_IXON = 'Andor iXon 888'

DEFAULT_PATH_IXON = '/home/fermikrb/luigi/imfit/testfiles/'
SENSOR_WIDTH_IXON = 1024

TPROBE_IXON = [40.0, 40.0]
ODSAT_IXON = [4.0, 4.0]
ISAT_FLUX_IXON = [275.0, 275.0]


########################################################################
######################     Fitting Parameters     ######################
########################################################################

DEFAULT_REGION = [[90, 125, 80, 80], 
                  [90, 125, 80, 80]
                 ]

FIT_FUNCTIONS = ['Gaussian', 'Bigaussian', 'Fermi-Dirac']


#######################################################################
########################    Imaging Parameters   ######################
#######################################################################

IMAGING_PATHS = ['Axial', 'Axial (Low Mag.)', 'Vertical']
PIXEL_SIZES = [1.8, 31.8, 1.29]
CSATEFF = [900.0, 1.0E6, 1400.0]

FRAMESEQUENCE = ['shadow', 'light', 'dark']
ATOM_NAMES = ['K','RB','KRB']
NATOMS = 2
WAVELENGTHS = [767.0E-9, 780.0E-9]


VERT_TRAP_ANGLE = 45.0 # Degrees

#######################################################################
########################  File System Parameters  #####################
#######################################################################

import datetime
now = datetime.datetime.now()

with open('./lib/ip.txt') as f:
	ip_str = f.read(100)
DEFAULT_PATH = now.strftime('//'+ip_str+'/krbdata/data/%Y/%m/%Y%m%d/') # PolarKRB's IP address


DEFAULT_PATH_PI = DEFAULT_PATH + 'ximea/'
DEFAULT_PATH_IXON = DEFAULT_PATH + 'Andor/'
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

SENSOR_WIDTH_IXON = 1024

TPROBE_IXON = [40.0, 40.0]
ODSAT_IXON = [4.0, 4.0]
ISAT_FLUX_IXON = [275.0, 275.0]

#######################################################################
########################        Ximea xiQ        ######################
#######################################################################

CAMERA_NAME_XIMEA = 'Ximea xiQ'


########################################################################
######################     Fitting Parameters     ######################
########################################################################

DEFAULT_REGION_XIMEA = [[200, 400, 300, 300], 
                  [200, 125, 300, 300]]
DEFAULT_REGION_IXON = [[130, 30, 100, 100], 
                  [130, 30, 100, 100]]
DEFAULT_REGION = {
	'XIMEA' : DEFAULT_REGION_XIMEA,
	'IXON' : DEFAULT_REGION_IXON
}

FIT_FUNCTIONS = ['Gaussian', 'Rotated Gaussian', 'Twisted Gaussian', 'Bigaussian', 'Fermi-Dirac', 'Vertical BandMap']
WORKSHEET_NAMES = ['Gauss1', 'Gauss1', 'Gauss1', 'Gauss2', 'FermiDirac', 'BandMapV']

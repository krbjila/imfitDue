import datetime

# from matplotlib.pyplot import autoscale
now = datetime.datetime.now()

DEFAULT_PATH = now.strftime(
    "K:/data/%Y/%m/%Y%m%d/"
)  # The dataserver is mounted to K: on the Windows computers.
FILESEP = "/"

AUTOSCALE_MIN = 2  # percentile
AUTOSCALE_HEADROOM = 1.1  # factor above max

# FIT_FUNCTIONS = ['Gaussian w/ Gradient', 'Gaussian', 'Rotated Gaussian', 'Twisted Gaussian', 'Bigaussian', 'Fermi-Dirac', 'Vertical BandMap']
FIT_FUNCTIONS = [
    "Gaussian w/ Gradient",
    "Gaussian",
    "Rotated Gaussian",
    "Twisted Gaussian",
    "Gaussian Fixed",
    "Bigaussian",
    "Fermi-Dirac",
    "Fermi-Dirac 2D",
    "Fermi-Dirac 2D Int",
    "Thomas-Fermi",
    "Integrate",
]
KRB_FIT_FUNCTIONS = [
    "Gaussian w/ Gradient",
    "Gaussian",
    "Rotated Gaussian",
    "Twisted Gaussian",
    "Gaussian Fixed",
    "Fermi-Dirac",
    "Fermi-Dirac 2D",
    "Fermi-Dirac 2D Int",
    "Integrate",
]
NONTWISTED_FIT_FUNCTIONS = [
    "Gaussian w/ Gradient",
    "Gaussian",
    "Bigaussian",
    "Fermi-Dirac",
    "Fermi-Dirac 2D Int",
    "Thomas-Fermi",
    "Integrate",
]
NONTWISTED_KRB_FIT_FUNCTIONS = [
    "Gaussian w/ Gradient",
    "Gaussian",
    "Fermi-Dirac",
    "Fermi-Dirac 2D Int",
    "Integrate",
]
WORKSHEET_NAMES = [
    "GaussGrad",
    "GaussGrad",
    "GaussGrad",
    "GaussGrad",
    "GaussGrad",
    "Gauss2",
    "FermiDirac",
    "FD2D",
    "FD2D",
    "ThomasFermi",
    "Integrated",
]

DEFAULT_MODE = "Side iXon"

MAX_OD_FIT = 1e9

CSAT = {
    "axial": {"K": 2970, "Rb": 2882},  # not calibrated for a while
    "side": {"K": 2955, "Rb": 3276},  # {"K": 2072, "Rb": 2201} 12/10/2024; {"K": 2955, "Rb": 3276} calibrated 12/02/2024
    "vertical": {"K": 490, "Rb": 1300},  # Calibrated 12/17/24
}

# TODO: Check these!
# Numerical aperture
NA = {"axial": 0.12, "side": 0.20, "vertical": 0.5}

TOP_CAMERA_ANGLE = 27.7 #27.7 degree, calibrated on 12.12.2024

# Pixel size (um)
PX_SIZE = {
    "axial": 2.58,
    "side": 1.718,
    "vertical": 0.92,
}  # vertical pixel size calibrated on Dec.01.2024

# TODO: Check these!
# Resonant cross section at I/Isat = 0, in um^2
SIGMA_0_K = 0.5 * 0.2807  # From Tiecke 40K data, assuming pi polarization
SIGMA_0_Rb = 0.5 * 0.2907  # From Steck 87Rb data, table 7, assuming pi polarization
SIGMA_0 = {
    "K": SIGMA_0_K,
    "Rb": SIGMA_0_Rb,
    "|0,0>": SIGMA_0_K,
    "|1,0>": SIGMA_0_K,
}

# Transfer efficiency
EFF_K = 1
EFF_Rb = 1
EFF_GSM = 0.86 * 0.79   # Feshbach imaging x 4.5 kV/cm STIRAP eff
#EFF_GSM = 0.83 * 0.9 # Feshbach imaging x at zero field
EFF = {
    "K": EFF_K,
    "Rb": EFF_Rb,
    "|0,0>": EFF_GSM,
    "|1,0>": EFF_GSM,
}

from collections import OrderedDict

IMFIT_MODES = OrderedDict(
    [
        (
            "Side iXon",
            {
                ### Mode 1 - Side iXon
                "Default Path": DEFAULT_PATH + "Andor/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["side"],
                "Species": ["K", "Rb"],
                "Image Path": "Side",
                "Default Region": [[94, 242, 200, 200], [94, 396, 150, 150]],
                "Extension Filter": "*.npz",
                "Fit Functions": NONTWISTED_FIT_FUNCTIONS,
                "Enforce same fit for both": False,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "K": {
                        "Shadow": 0,
                        "Light": 1,
                        "Dark": 2,
                    },
                    "Rb": {"Shadow": 3, "Light": 4, "Dark": 5},
                },
                "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
                "CSat": CSAT["side"],  # Unbinned effective C_sat
                "NA": NA["side"],
            },
        ),
        (
            "Side iXon SG", # Stern-Gerlach imaging - slightly hacky imaging of other K location on Rb frame.
            {
                ### Mode 1 - Side iXon
                "Default Path": DEFAULT_PATH + "Andor/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["side"],
                "Species": ["K", "K"],
                "Image Path": "Side",
                "Default Region": [[105, 177, 200, 50], [105, 225, 200, 50]],
                "Extension Filter": "*.npz",
                "Fit Functions": NONTWISTED_FIT_FUNCTIONS,
                "Enforce same fit for both": False,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "K": {
                        "Shadow": 0,
                        "Light": 3,
                        "Dark": 5,
                    },
                    "K": {"Shadow": 0, "Light": 3, "Dark": 5},
                },
                "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
                "CSat": {"K": 1200, "Rb": 1200},#CSAT["side"],  # Unbinned effective C_sat
                "NA": NA["side"],
            },
        ),
        (
            "Axial iXon",
            {
                ### Mode 1 - Axial iXon
                "Default Path": DEFAULT_PATH + "Andor/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["axial"],
                "Species": ["K", "Rb"],
                "Image Path": "Axial",
                "Default Region": [[140, 215, 250, 250], [130, 355, 350, 350]],
                "Extension Filter": "*.npz",
                "Fit Functions": FIT_FUNCTIONS,
                "Enforce same fit for both": False,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "K": {
                        "Shadow": 0,
                        "Light": 1,
                        "Dark": 2,
                    },
                    "Rb": {"Shadow": 3, "Light": 4, "Dark": 5},
                },
                "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
                "CSat": CSAT["axial"],  # Unbinned effective C_sat
                "NA": NA["axial"],
            },
        ),
        (
            "Side iXon Molecules In Situ",
            {
                ### Mode 4 - Side iXon Molecules
                "Default Path": DEFAULT_PATH + "MoleculeInSituFK/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["side"],
                "Species": ["|0,0>", "|1,0>"],
                "Image Path": "Side",
                "Default Region": [[186, 376, 150, 50], [186, 376, 150, 50]],
                "Extension Filter": "*.npz",
                "Fit Functions": NONTWISTED_KRB_FIT_FUNCTIONS,
                "Enforce same fit for both": True,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "|0,0>": {
                        "Shadow": 0,
                        "Light": 1,
                        "Dark": 2,
                    },
                    "|1,0>": {"Shadow": 3, "Light": 4, "Dark": 5},
                },
                "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
                "CSat": {
                    "|0,0>": CSAT["side"]["K"],
                    "|1,0>": CSAT["side"]["K"],
                },  # Unbinned effective C_sat,
                "NA": NA["side"],
                "Allow fit both states": True,
            },
        ),
        (
            "Side iXon Molecules ToF",
            {
                ### Mode 4 - Side iXon Molecules
                "Default Path": DEFAULT_PATH + "MoleculeInSituFK/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["side"],
                "Species": ["|0,0>", "|1,0>"],
                "Image Path": "Side",
                "Default Region": [[186, 400, 200, 400], [186, 400, 50, 50]],
                "Extension Filter": "*.npz",
                "Fit Functions": NONTWISTED_KRB_FIT_FUNCTIONS,
                "Enforce same fit for both": True,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "|0,0>": {
                        "Shadow": 3,
                        "Light": 4,
                        "Dark": 5,
                    },
                    "|1,0>": {"Shadow": 0, "Light": 1, "Dark": 2},
                },
                "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
                "CSat": {
                    "|0,0>": CSAT["axial"]["K"],
                    "|1,0>": CSAT["axial"]["K"],
                },  # Unbinned effective C_sat,
                "NA": NA["side"],
                "Allow fit both states": False,
            },
        ),
        # (
        #     "Axial iXon Molecules ToF",
        #     {
        #         ### Mode 1.1 - Axial iXon
        #         "Default Path": DEFAULT_PATH + "KRbFK/",
        #         "Default Suffix": "krbfk_{}.npz",
        #         "Pixel Size": PX_SIZE["axial"],
        #         "Species": ["|0,0>", "|1,0>"],
        #         "Image Path": "Axial",
        #         "Default Region": [[135, 150, 100, 100], [135, 150, 100, 100]],
        #         "Extension Filter": "*.npz",
        #         "Fit Functions": NONTWISTED_FIT_FUNCTIONS,
        #         "Enforce same fit for both": True,
        #         "Auto Detect Binning": True,
        #         "Array Width": 512,
        #         "Number of Frames": 4,
        #         "Frame Order": {
        #             "|0,0>": {
        #                 "Shadow": 0,
        #                 "Light": 2,
        #                 "Dark": 1,
        #             },
        #             "|1,0>": {"Shadow": 0, "Light": 2, "Dark": 1},
        #         },
        #         "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
        #         "CSat": {
        #             "|0,0>": CSAT["axial"]["K"],
        #             "|1,0>": CSAT["axial"]["K"],
        #         },  # Unbinned effective C_sat,
        #         "NA": NA["axial"],
        #         "Allow fit both states": True,
        #     },
        # ),
        # (
        #     "Axial iXon Molecules In Situ",
        #     {
        #         ### Mode 3 - Axial iXon Molecules
        #         "Default Path": DEFAULT_PATH + "MoleculeInSituFK/",
        #         "Default Suffix": "ixon_{}.npz",
        #         "Pixel Size": PX_SIZE["axial"],
        #         "Species": ["|0,0>", "|1,0>"],
        #         "Image Path": "Axial",
        #         "Default Region": [[125, 348, 100, 75], [125, 348, 100, 75]],
        #         "Extension Filter": "*.npz",
        #         "Fit Functions": KRB_FIT_FUNCTIONS,
        #         "Enforce same fit for both": True,
        #         "Auto Detect Binning": True,
        #         "Array Width": 512,
        #         "Number of Frames": 6,
        #         "Frame Order": {
        #             "|0,0>": {
        #                 "Shadow": 0,
        #                 "Light": 1,
        #                 "Dark": 2,
        #             },
        #             "|1,0>": {"Shadow": 3, "Light": 4, "Dark": 5},
        #         },
        #         "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
        #         "CSat": {
        #             "|0,0>": CSAT["axial"]["K"],
        #             "|1,0>": CSAT["axial"]["K"],
        #         },  # Unbinned effective C_sat,
        #         "NA": NA["axial"],
        #     },
        # ),
        # ('Axial iXon Molecules In Situ (CSV)', {
        #     ### Mode 3 - Axial iXon Molecules
        #     'Default Path': DEFAULT_PATH + 'MoleculeInSituFK/',
        #     'Default Suffix': 'ixon_{}.csv',
        #     'Pixel Size': PX_SIZE["axial"],
        #     'Species': ['|0,0>', '|1,0>'],
        #     'Image Path': 'Axial',
        #     'Default Region': [[272, 302, 60, 25],
        #               [272, 302, 60, 25]],
        #     'Extension Filter': '*.csv',
        #     'Fit Functions': KRB_FIT_FUNCTIONS,
        #     'Enforce same fit for both': True,
        #             'Auto Detect Binning': True,
        #     'Array Width': 512,
        #     'Number of Frames': 6,
        #     'Frame Order': {
        #         '|0,0>': {
        #             'Shadow': 0,
        #             'Light': 1,
        #             'Dark': 2,
        #         },
        #         '|1,0>': {
        #             'Shadow': 3,
        #             'Light': 4,
        #             'Dark': 5
        #         }
        #     },
        #     'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        #     'CSat': {'|0,0>': CSAT["axial"]["K"], '|1,0>': CSAT["axial"]["K"]}, # Unbinned effective C_sat,
        #     'NA': NA["axial"],
        # }),
        # (
        #     "Side iXon Molecules 4 Frame",
        #     {
        #         ### Mode 4 - Side iXon Molecules
        #         "Default Path": DEFAULT_PATH + "MoleculeFK4/",
        #         "Default Suffix": "ixon_{}.npz",
        #         "Pixel Size": PX_SIZE["side"],
        #         "Species": ["|0,0>", "|1,0>"],
        #         "Image Path": "Axial",
        #         "Default Region": [[115, 370, 80, 50], [115, 368, 80, 50]],
        #         "Extension Filter": "*.npz",
        #         "Fit Functions": KRB_FIT_FUNCTIONS,
        #         "Enforce same fit for both": True,
        #         "Auto Detect Binning": True,
        #         "Array Width": 256,
        #         "Number of Frames": 8,
        #         "Frame Order": {
        #             "|0,0>": {
        #                 "Shadow": 2,
        #                 "Light": 6,
        #                 "Dark": 3,
        #             },
        #             "|1,0>": {"Shadow": 4, "Light": 6, "Dark": 5},
        #         },
        #         "Fit angle": 0.0,  # Deg, Twisted Gaussian fit
        #         "CSat": {
        #             "|0,0>": CSAT["side"]["K"],
        #             "|1,0>": CSAT["side"]["K"],
        #         },  # Unbinned effective C_sat,
        #         "NA": NA["side"],
        #     },
        # ),
        # ('Pixelfly Test', {
        #     ### PLACEHOLDER
        #     ### TODO: put in finalized values
        #     'Default Path': DEFAULT_PATH + 'Pixelfly/',
        #     'Default Suffix': 'pixelfly_{}.npz',
        #     'Pixel Size': 2.58,
        #     'Species': ['K', 'Rb'],
        #     'Image Path': 'Axial',
        #     'Default Region': [[150, 220, 250, 250],
        #               [150, 350, 300, 300]],
        #     'Extension Filter': '*.npz',
        #     'Fit Functions': FIT_FUNCTIONS,
        #     'Enforce same fit for both': False,
        #     'Auto Detect Binning': False,
        #     'Number of Frames': 6,
        #     'Frame Order': {
        #         'K': {
        #             'Shadow': 0,
        #             'Light': 2,
        #             'Dark': 4,
        #         },
        #         'Rb': {
        #             'Shadow': 1,
        #             'Light': 3,
        #             'Dark': 5
        #         }
        #     },
        #     'Fit angle': 0.0, # Deg, Twisted Gaussian fit
        #     'CSat': {'K': 19e3 / 4.0, 'Rb': 19e3 / 4.0}, # Unbinned effective C_sat,
        #     'NA': NA["side"],
        # }),
        (
            "Vertical iXon",
            {
                ### Mode 2 - Vertical iXon
                "Default Path": DEFAULT_PATH + "Andor_Vertical/",
                "Default Suffix": "twospecies_{}.npz",
                "Pixel Size": PX_SIZE["vertical"],
                "Species": ["K", "Rb"],
                "Image Path": "Vertical",
                "Default Region": [[128, 295, 150, 150], [128, 295, 150, 150]],
                "Extension Filter": "*.npz",
                "Fit Functions": FIT_FUNCTIONS,
                "Enforce same fit for both": False,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "K": {
                        "Shadow": 0,
                        "Light": 1,
                        "Dark": 2,
                    },
                    "Rb": {"Shadow": 3, "Light": 4, "Dark": 5},
                },
                "Fit angle": TOP_CAMERA_ANGLE,  # Deg, Twisted Gaussian fit
                "CSat": CSAT["vertical"],  # Unbinned effective C_sat,
                "NA": NA["vertical"],
            },
        ),
        (
            "Vertical iXon Molecules",
            {
                ### Mode 2 - Vertical iXon
                "Default Path": DEFAULT_PATH + "KRb_Vertical/",
                "Default Suffix": "ixon_{}.npz",
                "Pixel Size": PX_SIZE["vertical"],
                "Species": ["|0,0>", "|1,0>"],
                "Image Path": "Vertical",
                "Default Region": [[224, 641, 200, 200], [224, 641, 200, 200]],
                "Extension Filter": "*.npz",
                "Fit Functions": KRB_FIT_FUNCTIONS,
                "Enforce same fit for both": True,
                "Auto Detect Binning": True,
                "Array Width": 512,
                "Number of Frames": 6,
                "Frame Order": {
                    "|0,0>": {
                        "Shadow": 0,
                        "Light": 1,
                        "Dark": 2,
                    },
                    "|1,0>": {"Shadow": 3, "Light": 4, "Dark": 5},
                },
                "Fit angle": TOP_CAMERA_ANGLE,  # Deg, Twisted Gaussian fit
                "CSat": {
                    "|0,0>": CSAT["vertical"]["K"],
                    "|1,0>": CSAT["vertical"]["K"],
                },  # Unbinned effective C_sat,
                "NA": NA["vertical"],
                "Allow fit both states": True,
            },
        ),
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
    ]
)

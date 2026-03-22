import numpy as np
from lib.polylog import dilog
import mpmath as mp
from scipy.special import gamma
import lib.polylog as polylog_lib
from lib.imfitDefaults import FREQS, PX_SIZE, NAT_CONSTANTS


def azimGauss(p, x, y):
    ### Parameters: [offset, amplitude, wx, dODdx]
    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -(x ** 2.0) / (2.0 * p[2] ** 2.0)
        )
        + p[3] * x
        - y
    )


def gaussian(p, r, y, mask_above=np.inf):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    mask = np.where(y < mask_above, 1, 0)

    XR = X * np.cos(p[6]) - Y * np.sin(p[6])
    YR = X * np.sin(p[6]) + Y * np.cos(p[6])

    x0R = p[2] * np.cos(p[6]) - p[4] * np.sin(p[6])
    y0R = p[2] * np.sin(p[6]) + p[4] * np.cos(p[6])

    return (
        np.ravel(
            p[0]
            + p[1]
            * np.exp(
                -((XR - x0R) ** 2.0) / (2.0 * p[3] ** 2.0)
                - (YR - y0R) ** 2.0 / (2.0 * p[5] ** 2.0)
            )
            - y
        )
        * mask.ravel()
    )


# Gaussian with Gradient; Mask above a certain radius
def gaussian_mask_sigma(
    p,
    r,
    y,
    mask_radius_x=1e-6,  # should be 1E-6 for default, 10 for testing
    mask_radius_y=1e-6,
    x0_mask=0,
    y0_mask=0,
    theta_mask=0,
):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    XR = X * np.cos(p[6]) - Y * np.sin(p[6])
    YR = X * np.sin(p[6]) + Y * np.cos(p[6])

    x0R = p[2] * np.cos(p[6]) - p[4] * np.sin(p[6])
    y0R = p[2] * np.sin(p[6]) + p[4] * np.cos(p[6])

    z = np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((XR - x0R) ** 2.0) / (2.0 * p[3] ** 2.0)
            - (YR - y0R) ** 2.0 / (2.0 * p[5] ** 2.0)
        )
        + p[7] * (XR - x0R)
        + p[8] * (YR - y0R)
    )

    XR_mask = (X - x0_mask) * np.cos(theta_mask) - (Y - y0_mask) * np.sin(theta_mask)
    YR_mask = (X - x0_mask) * np.sin(theta_mask) + (Y - y0_mask) * np.cos(theta_mask)

    QR = (XR_mask) ** 2 / mask_radius_x**2 + (YR_mask) ** 2 / mask_radius_y**2

    mask = np.where(QR > 1, 1, 0)

    return (z - np.ravel(y)) * mask.ravel()


# Gaussian with Gradient; Mask above a certain radius
def gaussian_mask_sigma_no_rot(
    p,
    r,
    y,
    mask_radius_x=1e-6,  # should be 1E-6 for default, 10 for testing
    mask_radius_y=1e-6,
    x0_mask=0,
    y0_mask=0,
):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    z = np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((X - p[2]) ** 2.0) / (2.0 * p[3] ** 2.0)
            - ((Y - p[4]) ** 2.0) / (2.0 * p[5] ** 2.0)
        )
        + p[6] * (X - p[2])
        + p[7] * (Y - p[4])
    )

    X_mask = X - x0_mask
    Y_mask = Y - y0_mask

    QR = (X_mask) ** 2 / mask_radius_x**2 + (Y_mask) ** 2 / mask_radius_y**2

    mask = np.where(QR > 1, 1, 0)

    return (z - np.ravel(y)) * mask.ravel()

# Gaussian with Gradient; Mask above a certain radius
def gaussian_mask_sigma_1D(
    p,
    r,
    y,
    mask_radius_x=1e-6,  # should be 1E-6 for default
    x0_mask=0,
):
    ### Parameters: [offset, amplitude, x0, wx, dODdx]
    X = r

    z = np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((X - p[2]) ** 2.0) / (2.0 * p[3] ** 2.0)
        )
        + p[4] * (X - p[2])
    )

    X_mask = X - x0_mask

    QR = (X_mask) ** 2 / mask_radius_x**2
    mask = np.where(QR > 1, 1, 0)

    return (z - np.ravel(y)) * mask.ravel()


def gaussianGradient(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta, dODdx, dODdy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    XR = X * np.cos(p[6]) - Y * np.sin(p[6])
    YR = X * np.sin(p[6]) + Y * np.cos(p[6])

    x0R = p[2] * np.cos(p[6]) - p[4] * np.sin(p[6])
    y0R = p[2] * np.sin(p[6]) + p[4] * np.cos(p[6])

    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((XR - x0R) ** 2.0) / (2.0 * p[3] ** 2.0)
            - (YR - y0R) ** 2.0 / (2.0 * p[5] ** 2.0)
        )
        + p[7] * (XR - x0R)
        + p[8] * (YR - y0R)
        - y
    )


def gaussianNoRot(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((X - p[2]) ** 2.0) / (2.0 * p[3] ** 2.0)
            - (Y - p[4]) ** 2.0 / (2.0 * p[5] ** 2.0)
        )
        - y
    )


def gaussianNoRotGradient(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, dODdx, dODdy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((X - p[2]) ** 2.0) / (2.0 * p[3] ** 2.0)
            - (Y - p[4]) ** 2.0 / (2.0 * p[5] ** 2.0)
        )
        + p[6] * X
        + p[7] * Y
        - y
    )


def gaussianNoRotTwist(p, r, y, angle):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)
    theta = angle * np.pi / 180.0
    XX = X * np.cos(theta) - Y * np.sin(theta)
    YY = Y * np.cos(theta) + X * np.sin(theta)

    x0R = p[2] * np.cos(theta) - p[4] * np.sin(theta)
    y0R = p[2] * np.sin(theta) + p[4] * np.cos(theta)

    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((XX - x0R) ** 2.0) / (2.0 * p[3] ** 2.0)
            - (YY - y0R) ** 2.0 / (2.0 * p[5] ** 2.0)
        )
        - y
    )


def doubleGaussian(p, r, y):
    ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    return np.ravel(
        p[0]
        + p[1]
        * np.exp(
            -((X - p[7]) ** 2.0) / (2 * p[2] ** 2.0)
            - (Y - p[8]) ** 2.0 / (2 * p[3] ** 2.0)
        )
        + p[4]
        * np.exp(
            -((X - p[7]) ** 2.0) / (2 * p[5] ** 2.0)
            - (Y - p[8]) ** 2.0 / (2 * p[6] ** 2.0)
        )
        - y
    )


def thomasFermi(p, r, y):
    ### Parameters: [offset, ampTF, x0, rx, y0, ry, ampGauss, wx, wy]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    if isinstance(y, int):
        y = np.zeros(X.shape)

    return np.ravel(
        p[0]
        + p[1]
        # 3/2 power for integrated profile
        * np.power(
            np.maximum(
                1 - (X - p[2]) ** 2.0 / p[3] ** 2.0 - (Y - p[4]) ** 2.0 / p[5] ** 2.0, 0
            ),
            3 / 2,
        )
        + p[6]
        * np.exp(
            -(
                (X - p[2]) ** 2.0 / (2 * p[7] ** 2.0)
                + (Y - p[4]) ** 2.0 / (2 * p[8] ** 2.0)
            )
        )
        - y
    )


def fermiDirac(p, r, y, mask_above=np.inf):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]

    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    if isinstance(y, int):
        y = np.zeros(X.shape)

    mask = np.where(y < mask_above, 1, 0)

    return (
        p[0]
        - p[1]
        * dilog(
            -np.exp(
                p[6]
                - (X.ravel() - p[2]) ** 2.0 / (2 * p[3] ** 2.0)
                - (Y.ravel() - p[4]) ** 2 / (2 * p[5] ** 2)
            )
        )
        - y.ravel()
    ) * mask.ravel()

def fermiDirac_fixed_bemu(p, r, y, mask_above=np.inf, N0 = 1E3, TOF = 0,
                          omega_x = 2 * np.pi * FREQS["fx"],
                          omega_y = 2 * np.pi * FREQS["fy"],
                          omega_z = 2 * np.pi * FREQS["fz"],
                          pxsz_um = PX_SIZE["side"] * 1E-6,
                          mass = 40 * NAT_CONSTANTS["amu2kg"]):
    """
    Fermi-Dirac distribution with fixed chemical potential (betamu) to match a given number of particles N0.
    This is a 3D version that takes into account the TOF and the trap frequencies.
    """
    ### Parameters: [offset, amplitude, x0, wx, y0, wy]
    # Take N as a constraint
    ### constants
    hbar = NAT_CONSTANTS["hbar"] # J s
    kB = NAT_CONSTANTS["kB"]  # J/K

    omega_bar = (omega_x * omega_y * omega_z) ** (1 / 3)

    Tx = mass * omega_x**2 * (p[3] * pxsz_um)**2 / (1 + omega_x**2 * TOF**2) / kB
    Ty = mass * omega_y**2 * (p[5] * pxsz_um)**2 / (1 + omega_y**2 * TOF**2) / kB
    T_avg = (Tx**2 * Ty) ** (1/3) # take geometric mean of Tx and Ty
    
    # Find where betamu for the given T gives us the right number of particles
    betamu_range = np.linspace(-10, 20, 5000)
    N_checker_3D = (
                    (kB * T_avg / (hbar * omega_bar)) ** 3
                    ) * polylog_lib.fermi_poly3(betamu_range)
    
    N_diff_3D = np.abs(N_checker_3D - N0)
    ind_N_3D = np.where(N_diff_3D == np.min(N_diff_3D))[0]
    betamu_3D = betamu_range[ind_N_3D[0]]

    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    if isinstance(y, int):
        y = np.zeros(X.shape)

    mask = np.where(y < mask_above, 1, 0)

    return (
        p[0]
        - p[1]
        * dilog(
            -np.exp(
                betamu_3D
                - (X.ravel() - p[2]) ** 2.0 / (2 * p[3] ** 2.0)
                - (Y.ravel() - p[4]) ** 2 / (2 * p[5] ** 2)
            )
        )
        - y.ravel()
    ) * mask.ravel()


def fermiDirac2D(p, r, y, angle):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]

    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)
    theta = angle * np.pi / 180.0
    XX = X * np.cos(theta) - Y * np.sin(theta)
    YY = Y * np.cos(theta) + X * np.sin(theta)

    x0R = p[2] * np.cos(theta) - p[4] * np.sin(theta)
    y0R = p[2] * np.sin(theta) + p[4] * np.cos(theta)

    if isinstance(y, int) or isinstance(y, float):
        y = np.zeros(X.shape)

    return (
        p[0]
        + p[1]
        * np.log(
            1
            + np.exp(
                p[6]
                - (XX.ravel() - x0R) ** 2.0 / (2 * p[3] ** 2.0)
                - (YY.ravel() - y0R) ** 2.0 / (2 * p[5] ** 2.0)
            )
        )
        - y.ravel()
    )


def gaussian1D(p, r, y):
    ### Parameters: [offset, amplitude, x0, sigma, gradient]
    return p[0] + p[1] * np.exp(-((r - p[2]) ** 2.0) / (2 * p[3] ** 2.0)) + p[4] * r - y


polylog = mp.fp.polylog
polylog_np = np.vectorize(lambda n, x: np.real(mp.fp.polylog(n, x)))


def limexp(n, x):
    return np.where(x < 30, np.real(polylog_np(n, -np.exp(x))), -(x**n) / gamma(n + 1))


def fermiDirac2Dint(p, r, y):
    ### Parameters: [offset, amplitude, x0, sigma, q, gradient]

    return (
        p[0]
        - p[1] * limexp(3 / 2, p[4] - (r - p[2]) ** 2 / (2 * p[3] ** 2))
        + p[5] * r
        - y
    )


def bandmapV(p, r, y, imageDetails):
    from scipy.special import erf as erf

    mList = [40.0, 87.0, 127.0]

    atom = imageDetails[0]
    pxl = imageDetails[2]

    ################ Fixed parameters
    k = 2.0 * np.pi / (1.064e-6)
    hbar = 1.0545718e-34
    u = 1.66053904e-27
    m = mList[atom] * u
    TOF = imageDetails[1] * 1e-3

    Delta = 2.0 * hbar * k / m * TOF * 1e6 / pxl  ### Delta is in pixels

    ### Parameters: [offset, A0, A1, A2, wy, yc, wx, xc]

    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    B1 = (
        p[1]
        / 2.0
        * (
            erf((Delta + 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((Delta - 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
        )
    )
    B2 = (
        p[2]
        / 2.0
        * (
            erf((-Delta + 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((-Delta - 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((2.0 * Delta + 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((2.0 * Delta - 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
        )
    )
    B3 = (
        p[3]
        / 2.0
        * (
            erf((-2.0 * Delta + 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((-2.0 * Delta - 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((3.0 * Delta + 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
            + erf((3.0 * Delta - 2.0 * (Y - p[5])) / (2.0 * np.sqrt(2.0) * p[4]))
        )
    )

    GX = np.exp(-((X - p[7]) ** 2.0) / (2.0 * p[6] ** 2.0))

    return np.ravel((B1 + B2 + B3) * GX + p[0] - y)

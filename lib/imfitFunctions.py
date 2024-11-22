import numpy as np
from lib.polylog import dilog
import mpmath as mp
from scipy.special import gamma


def gaussian(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta]
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
        - y
    )


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


def fermiDirac(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]

    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    if isinstance(y, int):
        y = np.zeros(X.shape)

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
    )


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

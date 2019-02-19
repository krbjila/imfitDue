import numpy as np
from polylog import dilog

def gaussian(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, theta]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis, yaxis)

    XR = X*np.cos(p[6]) - Y*np.sin(p[6])
    YR = X*np.sin(p[6]) + Y*np.cos(p[6])

    x0R = p[2]*np.cos(p[6]) - p[4]*np.sin(p[6])
    y0R = p[2]*np.sin(p[6]) + p[4]*np.cos(p[6])

    return np.ravel(p[0] + p[1]*np.exp( -(XR-x0R)**2.0/(2*p[3]**2.0)  -(YR-y0R)**2.0/(2*p[5]**2.0)  ) - y)


def doubleGaussian(p,r,y):
    ### Parameters: [offset, Amp1, wx1, wy1, Amp2, wx2, wy2, x0, y0]
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis,yaxis)

    return np.ravel(p[0] + p[1]*np.exp(-(X-p[7])**2.0/(2*p[2]**2.0) -(Y-p[8])**2.0/(2*p[3]**2.0))\
            + p[4]*np.exp(-(X-p[7])**2.0/(2*p[5]**2.0) -(Y-p[8])**2.0/(2*p[6]**2.0)) - y)

def fermiDirac(p, r, y):
    ### Parameters: [offset, amplitude, x0, wx, y0, wy, q]
    
    xaxis = r[0]
    yaxis = r[1]

    X, Y = np.meshgrid(xaxis,yaxis)

    
    if isinstance(y,int):
        y = np.zeros(X.shape)

    return p[0] - p[1]*dilog(-np.exp(p[6]-(X.ravel()-p[2])**2.0/(2*p[3]**2.0)-(Y.ravel()-p[4])**2/(2*p[5]**2)))-y.ravel()
    


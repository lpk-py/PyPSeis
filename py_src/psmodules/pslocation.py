'''

# A python module for microseismic location using
# arrival-time picking based method, and
# migration style method, i.e semblance, linear stacking, semblance weighted stacking  
#
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

import os
from math import floor
import numpy as np
from numpy import array, ones, diff, square, cumsum, zeros, linspace, append, exp,\
    savetxt, loadtxt, argmin, abs
import matplotlib.pyplot as plt
import psarray
import pssynthetic
import psraytrace
import pswavelet
import psplot
import pspicker
import pspdf


def main():
    # Generate square grid array
    geox, geoy, geoz = psarray.gridarray(49, 6000, 6000)

    # geox = linspace(0, 6000, 31)
    # geoy = linspace(0, 0, 31)
    # geoz = geoy + 2.0
    nt = len(geox)

    # Define source coordinates
    sourcex = array([3000])
    sourcey = array([3000])
    sourcez = array([3000])

    # Define geological model
    zlayer = array([0, 80, 180, 300, 550, 930, 1300,
                    1650, 2000, 2300, 2680, 2920, 3100])
    zlayer_ss103h = array([0, 10, 30, 110, 210, 315, 395, 568, 821, 923, 1177, 1270, 1331, 1383, 1593, 1656, 1746,
                           2155, 2301, 2488, 2701, 2863, 2974, 3114, 3200])

    # Define velocity model
    # P wave velocity
    vp = array([1800, 2000, 2200, 2500, 2400, 2700,
                3150, 2950, 3440, 3750, 4000, 4350, 4600])
    vp_ss103h = array([880.00, 1200.00, 1732.00, 1732.00, 1732.00, 2122.00, 2193.00, 2550.00, 2522.00,
                       2709.00, 2793.00, 2588.00, 2596.00, 3047.00, 3170.00, 3061.00, 2537.00, 3573.00, 3655.00, 4205.00, 4345.00,
                       4486.00, 4629.00, 4460.00, 4843.00])
    vs_ss103h = array([400, 545, 787.00, 787.00, 870.00, 1060.00, 1100.00, 1272.00, 1256.00, 1350.00,
                       1470.00, 1362.00, 1366.00, 1603.00, 1668.00, 1611.00, 1270.00, 1985.00, 2030.00, 2336.00, 2508.00, 2590.00,
                       2673.00, 2575.00, 2796.00])

    #  S wave velocity based on Castagna's rule
    vs = (vp - 1360) / 1.16

    # 3D passive seismic raytracing
    print("3D passive seismic raytracing example is running[Waiting...]")
    dg = 10
    ptimes, pthetas = psraytrace.raytracing(
        vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sourcex, sourcey, sourcez, geox, geoy, geoz)

    print("3D passive seismic raytracing completed[OK]")
    # print(pthetas)
    # print("Trave times:")
    # print(ptimes)

    # Generate wavelet
    # oscillator wavelet
    # osciwlet, tw = pswavelet.oscillator(0.0005, 65, 0.2555, 3, 3, 1, 80, 50)
    osciwlet, tw = pswavelet.oscillator(0.002, 65, 1.022, 3, 3, 1, 80, 50)

    # Generate synthetic microseismogram
    # surface acquisition
    ns = 1500
    dt = 0.002
    # # downhole acquisition
    # ns = 600
    # dt = 0.0005

    data = zeros((ns, nt), dtype='float32')
    att = ones((nt, 1), dtype='float32')
    syndata = pssynthetic.genSynNoiseFree(
        ns, nt, osciwlet, pthetas, ptimes, dt, att)

    # nt = 62
    # Plot synthetic traces
    # psplot.hseisplot(syndata, ns, nt)
    # psplot.vseisplot(syndata, ns, nt)

    # Pick and plot first arrival time of microseismic events

    # pick arrival times
    pics = []
    for i in range(nt):
        tr = syndata[:, i]
        pic = pspicker.merpicker(tr, ns, 20, 600, "False")
        pics.append(pic)
    pickers = array(pics, dtype='float32')
    pickers.shape = (len(pickers), 1)

    # plot pickers
    # psplot.hseispickplot(syndata, pickers, ns, nt)
    # psplot.vseispickplot(syndata, pickers, ns, nt)

    # Locate microseismic event
    x1 = 0
    x2 = 6000
    dx = 600
    nx = floor((x2 - x1) / dx)

    y1 = 0
    y2 = 6000
    dy = 600
    ny = floor((y2 - y1) / dy)

    z1 = 2900
    z2 = 3100
    dz = 20
    nz = floor((z2 - z1) / dz)

    n = 1
    t0 = 0
    sigma = 1

    minErrs = []
    for i in range(x1, x2, dx):
        for j in range(y1, y2, dy):
            sx = array([i])
            sy = array([j])
            sz = array([3000])
            tps, tetas = psraytrace.raytracing(
                vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sx, sy, sz, geox, geoy, geoz)
            tps = tps / dt

            minErr = timediff(tps, pickers)
            minErrs.append(minErr)

    # pdfs = array(pdfs)
    minErrs = array(minErrs)
    minErrIndex = argmin(minErrs)
    print(minErrIndex)
    ind = findIndex(nx, ny, 1, minErrIndex)

    print(ind)


def timediff(tp, tmp):
    """
    timediff(): microseismic event location method based on arrival time differences. 
    This method needs to pick first arrival times of microseismic event and 
    generally aims to process high signal-to-noise ratio.

    """
    tpdiff = abs(diff(tp, axis=0))
    tmpdiff = abs(diff(tmp, axis=0))

    temp = square(tpdiff - tmpdiff)
    sumErrs = cumsum(temp)
    minErr = sumErrs[len(sumErrs) - 1]
    return minErr


def time():
    """
    time(): microseismic event location method based on arrival times. 
    This method needs to pick first arrival times of microseismic event and 
    generally aims to process high signal-to-noise ratio.

    """
    print('time()')


def sws():
    """
    sws(): semblance weighted stacking for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def semblance():
    """
    linearstack(): semblance for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def linearstack():
    """
    linearstack(): linear stacking for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def findIndex(nx, ny, nz, ind):
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                value = k * nx * ny + i * ny + j
                if value == ind:
                    return np.array([i, j, k])

# This will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

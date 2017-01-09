'''

# A python module for wavelet generation
#
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

from numpy import square, linspace, exp, sin, cos, concatenate
from math import ceil, pi


def main():
    """
    Run an example :
    ricker(0.002, 65, 1.022)
    """
    import matplotlib.pyplot as plt

    # # rikcer wavelet
    # wlet, tw = ricker(0.002, 65, 1.022)
    # plt.plot(wlet)
    # plt.show()

    # oscillator wavelet
    wlet, tw = oscillator(0.002, 65, 1.022, 3, 3, 1, 80, 50)
    plt.plot(wlet)
    plt.show()


def oscillator(dt=0.002, fdom=65, tlength=1.022, rlratio=3, A1=3, A2=1, k1=80, k2=50):
    """

    Both the P and S wavelets have the form
    s(t) = A0 * sin(2πft) * exp(–kt)
    This is the response of a damped harmonic oscillator. Specific
    values for the parameters f and k are arbitrary; we have chosen
    f = 300 Hz, k = 80/s for the P arrivals, and f = 200 Hz, k = 50/s
    for the S arrivals.

    dt = desired temporal sample rate in seconds
    fdom = dominant frequency in Hz (default: 65 Hz)
    tlength = wavelet lenghth in seconds (default: 511*dt, i.e. power of 2)

    Zhengguang Zhao, 2016, 
    Modified from G.F. Margrave, May 1991

    """

    # create a time vector
    n = round(tlength / dt) + 1
    nzero = ceil((n + 1) / rlratio)  # zero time sample is here
    nr = n - nzero  # number of samples to the right of nzero
    nl = n - nr - 1  # number of samples to the left of nzero

    ltw = dt * linspace(-nl, 0, nl + 1)
    ltw.shape = (len(ltw), 1)
    rtw = dt * linspace(0, nr, nr)
    rtw.shape = (len(rtw), 1)

    tw = concatenate((ltw, rtw))

    # create wavelet
    lwavelet = -1 * A1 * sin(2 * pi * fdom * ltw) * exp(-1 * k1 * ltw)
    rwavelet = 1 * A2 * sin(2 * pi * fdom * rtw) * exp(-1 * k2 * rtw)

    wavelet = concatenate((lwavelet, rwavelet))

    # (1 - 2 * pf * square(tw)) * exp(-pf * square(tw))

    # normalize
    # generate a reference sinusoid at the dominant frequency
    wavelet = wavenorm(wavelet, tw, 1)

    return wavelet, tw


def ricker(dt=0.002, fdom=65, tlength=1.022):
    """

    dt = desired temporal sample rate in seconds
    fdom = dominant frequency in Hz (default: 65 Hz)
    tlength = wavelet lenghth in seconds (default: 511*dt, i.e. power of 2)

    Zhengguang Zhao, 2016, 
    Modified from G.F. Margrave, May 1991

    """

    # create a time vector
    n = round(tlength / dt) + 1
    nzero = ceil((n + 1) / 2)  # zero time sample is here
    nr = n - nzero  # number of samples to the right of nzero
    nl = n - nr - 1  # number of samples to the left of nzero

    tw = dt * linspace(-nl, nr, nr + nl + 1)
    tw.shape = (len(tw), 1)

    # create wavelet
    pf = pi * fdom
    pf = pf * pf

    wavelet = (1 - 2 * pf * square(tw)) * exp(-pf * square(tw))

    # normalize
    # generate a reference sinusoid at the dominant frequency
    wavelet = wavenorm(wavelet, tw, 1)

    return wavelet, tw


def wavenorm(w, tw, flag):
    """
    wnorm = wavenorm(w,tw,flag)

    WAVENORM normalizes a wavelet by one of several criteria.
    The choices are: (1) normalize the maximum absolut value to unity
    (2) normalize such that a sine wave at the dominant frequency passes
    with unit amplitude (3) normalize the rms amplitude to unity.

    w ... input wavelet
    tw ... time coordinate vector for w
    flag ... (1) normalize the maximum absolute value to unity
            (2) normalize such that a sine wave at the dominant frequency 
                passes with unit amplitude
            (3) normalize the rms amplitude to unity

    Zhengguang Zhao, 2016, 
    Modified from G.F. Margrave, May 1991

    """

    if flag == 1:
        wnorm = w / max(abs(w))

    return wnorm


# this will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

'''

# A python module for synthetic seismogram generation
#
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

from math import cos, pi, floor
from numpy import zeros, array, append


def main():
    """
    Example for generating one seismic trace

    """

    import matplotlib.pyplot as plt
    import pswavelet

    ns = 2000
    nt = 1
    wlet, tw = pswavelet.oscillator(0.002, 65, 1.022, 3, 3, 1, 80, 50)
    dip = array([pi / 3])
    t = array([2])
    dt = array([0.002])
    att = array([1.0])
    syndata = genSynNoiseFree(ns, nt, wlet, dip, t, dt, att)

    plt.plot(syndata)
    plt.show()


def genSynNoiseFree(ns, nt, wavelet, dip, time, dt, att):
    """
    generate synthetic seismogram
    data = trace data of a certain number of samples
    wavelet = wavelet used for synthetic data
    dip = dip angle vector for each source-receiver pair
    time = travel time vector in seconds
    dt = sample rate in seconds,
    att = attenuation coefficient for each source-receiver pair based on distance

    """
    syn = zeros((ns, nt), dtype='float32')
    win = len(wavelet)

    for i in range(nt):
        nstart = floor(time[i, :] / dt)
        trace = 0.0 - cos(dip[i, :]) * wavelet * att[i]
        syn[nstart:nstart + win, i] += trace[:, 0]

    return syn


# This will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

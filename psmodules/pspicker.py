'''

# A python module for picking P- and/or S-wave first arrivals
# of microseismic event
# 
# 
# (C) Zhengguang Zhao, 2016

'''
# !/Users/Uqer/anaconda/bin/python

import numpy as np
from numpy import square, where, empty, arange

import matplotlib.pyplot as plt


def main():
    from obspy.segy.segy import readSEGY
    filename = '/Users/Uqer/Dropbox/SEGY/dr01.sgy'
    st = readSEGY(filename)
    tr = st.traces[0].data
    picker = merpicker(tr, 512, 0.25, 20, 600, "True")
    print(picker)


def erpicker(data, ns, win, threshold, plot='False'):
    """
    First arrival time picker based on energy ratio method
    """

    er = np.zeros((ns - win,), dtype='float32')
    for i in range(win, ns - win):
        ltemp = square(data[i - win + 1: i])
        le = sum(ltemp)
        rtemp = square(data[i: i + win])
        re = sum(rtemp)
        er[i] = re / le
    ermax = max(er)
    # print(ermax)
    if threshold < ermax:
        ind = where(er > threshold)
        pic = ind[0][0]
    else:
        ind = where(er == ermax)
        pic = ind[0][0]

    if plot == 'True':
        fig = plt.figure(figsize=(12, 8))

        xcor = arange(1, ns + 1, 1)
        ylen = len(xcor)

        absSegyTraceData = empty([ns, 1])

        for i in range(ylen):
            absSegyTraceData[i] = abs(tr[i])
        ymax = max(absSegyTraceData)
        TraceData = empty([ns, 1])
        for i in range(ylen):
            TraceData[i] = tr[i] / ymax

        ax1 = plt.subplot(2, 1, 1)
        ax1.plot(xcor, TraceData, color="black")
        ax1.axis([1, ns, -1, 1])
        plt.hold(True)
        ax1.plot(pic, 0, marker="o", color="red")
        ax2 = plt.subplot(2, 1, 2)
        ax2.plot(er, color="blue")
        plt.show()

    return pic


def merpicker(data, ns, win, threshold, plot='true'):
    """
    First arrival time picker based on modified energy ratio method
    """
    mer = np.zeros((ns - win,), dtype='float32')
    for i in range(win, ns - win):
        ltemp = square(data[i - win + 1: i])
        le = sum(ltemp)
        rtemp = square(data[i: i + win])
        re = sum(rtemp)
        mer[i] = pow((re / le) * abs(data[i]), 3)
    mermax = max(mer)
    # print(mermax)
    if threshold < mermax:
        ind = where(mer > threshold)
        pic = ind[0][0]
    else:
        ind = where(mer == mermax)
        pic = ind[0][0]

    if plot == 'true':
        fig = plt.figure(figsize=(12, 8))

        xcor = arange(1, ns + 1, 1)
        ylen = len(xcor)

        absSegyTraceData = empty([ns, 1])

        for i in range(ylen):
            absSegyTraceData[i] = abs(data[i])
        ymax = max(absSegyTraceData)
        TraceData = empty([ns, 1])
        for i in range(ylen):
            TraceData[i] = data[i] / ymax

        ax1 = plt.subplot(2, 1, 1)
        ax1.plot(xcor, TraceData, color="black")
        ax1.axis([1, ns, -1, 1])
        plt.hold(True)
        ax1.plot(pic, 0, marker="o", color="red")
        ax2 = plt.subplot(2, 1, 2)
        ax2.plot(mer, color="blue")
        plt.show()

    return pic


# This will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

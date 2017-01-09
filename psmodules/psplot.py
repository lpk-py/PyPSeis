
"""

A python module for plotting seismic traces in varieties of ways, i.e.


"""


#
# psplotseis : A Python module for plotting seismic traces
#
# (C) Zhengguang Zhao, 2016
#


# !/Users/Uqer/anaconda/bin/python


import numpy as np


def main():

    filename = '/Users/Uqer/Dropbox/SEGY/dr01.sgy'

    # use obspy module to read SEGY filename
    from obspy import read
    st = read(filename)
    # print(st.traces[0].data)
    # st.plot()
    ns = st.traces[0].stats.npts
    data = st.traces
    nt = len(st.traces)
    vdata = np.array(data).transpose()

    imageSegy(vdata)
    imageSegy(data)
    wiggleHAll(data, nt, ns)
    wiggleVAll(data, nt, ns)

    # # use pssegy module to read SEGY file
    # import pssegy
    # pssegy.verbose = 1
    # Data, SH, STH = pssegy.readSegy(filename)
    # imageSegy(Data)
    # pssegy.wiggleVComp(Data, SH, 1, 0.5, 'blue', 0.2)


def imageSegy(Data):
    """
    imageSegy(Data)
    Image segy Data
    """
    import pylab
    pylab.imshow(Data, cmap='seismic')
    pylab.title('Variable Density Display')
    pylab.grid(True)
    pylab.show()


def hseispickplot(Data, Picks, ns, nt, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Horizontal Wiggle Display'):
    """
    hseispickplot(Data, Picks, nt,ns ...) display seismic traces and event picks horizontally
    Data = seismic data matrix of dimension (ns, nt)
    Picks = event picks matrix of dimension(nt,1)

    """

    import matplotlib.pyplot as plt
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[:, i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[:, i])

        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        plt.plot(t, 1 + i + normtrace * scale,
                 color=ucolor, linewidth=lwidth)
        plt.plot(Picks[i, 0], 1 + i,
                 marker="o", color="red")
        plt.hold(True)
        # for a in range(len(trace)):
        #     if (trace[a] < 0):
        #         trace[a] = 0
        # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    plt.title(title)
    plt.axis([0, ns, 0, nt + 1])

    plt.grid(True)
    plt.show()


def vseispickplot(Data, Picks, ns, nt, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Vertical Wiggle Display'):
    """
    vseispickplot(Data, Picks, nt,ns ...) display seismic traces and event picks vertically
    Data = seismic data matrix of dimension (ns, nt)
    Picks = event picks matrix of dimension(nt,1)

    """

    import matplotlib.pyplot as plt
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[:, i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[:, i])
        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        plt.plot(1 + i + normtrace * scale,
                 t, color=ucolor, linewidth=lwidth)
        plt.plot(1 + i, Picks[i, 0], marker="o", color="red")
        for a in range(len(trace)):
            if (trace[a] < 0):
                trace[a] = 0
            # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    plt.title(title)
    plt.axis([0, int(nt / skipt) + 1, ns, 1])
    plt.grid(True)
    plt.show()


def hseisplot(Data, ns, nt, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Horizontal Wiggle Display'):
    """
    hseisplot(Data, nt,ns ...) display seismic traces horizontally
    Data = seismic data matrix of dimension (ns, nt)
    """

    import matplotlib.pyplot as plt
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[:, i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[:, i])

        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        plt.plot(t, 1 + i + normtrace * scale,
                 color=ucolor, linewidth=lwidth)
        plt.hold(True)
        # for a in range(len(trace)):
        #     if (trace[a] < 0):
        #         trace[a] = 0
        # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    plt.title(title)
    plt.axis([0, ns, 0, nt + 1])

    plt.grid(True)
    plt.show()


def vseisplot(Data, ns, nt, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Vertical Wiggle Display'):
    """
    vseisplot(Data, nt,ns ...) display seismic traces vertically
    Data = seismic data matrix of dimension (ns, nt)

    """

    import matplotlib.pyplot as plt
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[:, i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[:, i])
        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        plt.plot(1 + i + normtrace * scale,
                 t, color=ucolor, linewidth=lwidth)
        for a in range(len(trace)):
            if (trace[a] < 0):
                trace[a] = 0
            # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    plt.title(title)
    plt.axis([0, int(nt / skipt) + 1, ns, 1])
    plt.grid(True)
    plt.show()


def wiggleVAll(Data, nt, ns, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Vertical Wiggle Display'):
    """
    wiggle(Data,SH)
    """

    import pylab
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[:, i])
        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        pylab.plot(1 + i + normtrace * scale,
                   t, color=ucolor, linewidth=lwidth)
        for a in range(len(trace)):
            if (trace[a] < 0):
                trace[a] = 0
            # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    pylab.title(title)
    pylab.axis([0, nt + 1, ns, 1])

    pylab.grid(True)
    pylab.show()


def wiggleVComp(Data, nt, ns, component=1, scale=0.5, ucolor='blue', lwidth=1.0, title="Vertical Wiggle Display"):
    """
    wiggle(Data,SH)
    """

    import pylab
    from numpy import abs

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    if (component != 3):
        for i in range(component, nt, 3):
            #		trace=zeros(SH['ns']+2)
            #		dtrace=Data[:,i]
            #		trace[1:SH['ns']]=Data[:,i]
            #		trace[SH['ns']+1]=0
            trace = np.array(Data[i])
            if trace.max() > abs(trace.min()):
                maxval = trace.max()
            else:
                maxval = abs(trace.min())

            trace[0] = 0
            trace[ns - 1] = 0

            pylab.plot(1 + i / 3 + trace / maxval * scale,
                       t, color=ucolor, linewidth=lwidth)
            for a in range(len(trace)):
                if (trace[a] < 0):
                    trace[a] = 0
                # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)
    else:
        for i in range(0, nt, 1):
            #		trace=zeros(SH['ns']+2)
            #		dtrace=Data[:,i]
            #		trace[1:SH['ns']]=Data[:,i]
            #		trace[SH['ns']+1]=0
            trace = Data[:, i]
            if trace.max() > abs(trace.min()):
                maxval = trace.max()
            else:
                maxval = abs(trace.min())

            trace[0] = 0
            trace[ns - 1] = 0

            # Display each component with different color
            if (i % 3 == 0):
                ucolor = 'green'
            elif (i % 3 == 1):
                ucolor = 'blue'
            else:
                ucolor = 'red'

            pylab.plot(1 + i / 3 + trace / maxval * scale,
                       t, color=ucolor, linewidth=lwidth)
            for a in range(len(trace)):
                if (trace[a] < 0):
                    trace[a] = 0

    pylab.title(title)
    pylab.axis([0, nt / 3 + 1, ns, 1])
    pylab.grid(True)
    pylab.show()


def wiggleHAll(Data, nt, ns, skipt=1, scale=0.5, ucolor='blue', lwidth=1.0, title='Horizontal Wiggle Display'):
    """
    wiggle(Data,SH)
    """

    import matplotlib.pyplot as plt
    from numpy import abs, empty

    # find max of each trace to normalization
    norm = maxval = empty((nt, 1), dtype='float32')
    for i in range(0, nt, skipt):
        trace = np.array(Data[i])
        if trace.max() > abs(trace.min()):
            maxval[i] = trace.max()
        else:
            maxval[i] = abs(trace.min())
    norm = maxval / maxval.max()

    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = np.array(Data[i])

        trace[0] = 0
        trace[ns - 1] = 0
        normtrace = trace / maxval[i] * norm[i]
        plt.plot(t, 1 + i + normtrace * scale,
                 color=ucolor, linewidth=lwidth)
        plt.hold(True)
        # for a in range(len(trace)):
        #     if (trace[a] < 0):
        #         trace[a] = 0
        # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    plt.title(title)
    plt.axis([0, ns, 0, nt + 1])

    plt.grid(True)
    plt.show()

# this will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

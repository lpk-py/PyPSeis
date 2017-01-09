'''

# A python module for generating microseismic surface arrays, i.e.
# star-shaped array
# patch array
# grid array (dense, sparse, square, circular, etc.)
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python


    
import math
import numpy as np

# Global parameters
X_LENGTH = 6000
Y_LENGTH = 6000
GEO_NUM = 1000
X_INNER = 3000
Y_INNER = 3000

## Generating Grid Array
# Case 1: Square Grid
def gridarray(geonum = 1000, xlen = 6000, ylen = 6000):

    ngeo = int(math.sqrt(geonum))
    GEO_X_NUM = ngeo
    GEO_Y_NUM = ngeo
    # print('Number of geophones in XY directions are [%d,%d]' %(GEO_X_NUM, GEO_Y_NUM))
    X_OFFSET = xlen/(GEO_X_NUM-1)
    # print ('Offset in X direction is %d' %(X_OFFSET))
    Y_OFFSET = ylen/(GEO_Y_NUM-1)
    # print ('Offset in Y direction is %d' %(Y_OFFSET))
    # print('Geophone offsets in XY directions are [%d, %d]' %(X_OFFSET, Y_OFFSET))
    total_geo_num = GEO_X_NUM*GEO_Y_NUM
    # print ('Total number of geophones used in Square Array is: %d' %(total_geo_num))
    geo_x = np.zeros((total_geo_num,), dtype = np.float32)
    geo_y = np.zeros((total_geo_num,), dtype = np.float32)
    geo_z = np.zeros((total_geo_num,), dtype = np.float32)

    for i in range(GEO_Y_NUM):
        for j in range(GEO_X_NUM):
            geo_x[i * GEO_X_NUM + j] = j * X_OFFSET
            geo_y[i * GEO_X_NUM + j] = i * Y_OFFSET
            geo_z[i * GEO_X_NUM + j] = 2.0
            print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' %(i * GEO_X_NUM + j, geo_x[i * GEO_X_NUM + j], geo_y[i * GEO_X_NUM + j], geo_z[i * GEO_X_NUM + j]))

    # for i in range(len(geox)):
    #     print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' % (i , geox[i], geoy[i], geoz[i]))

    return geo_x, geo_y, geo_z

# Generating downhole array
def dharray(ngeo = 12, x= 0, y = 0, bottomz = 3000, voffset = 25):

    total_geo_num = ngeo
    # print ('Total number of geophones used in Square Array is: %d' %(total_geo_num))
    geo_x = np.zeros((total_geo_num,), dtype = np.float32)
    geo_y = np.zeros((total_geo_num,), dtype = np.float32)
    geo_z = np.zeros((total_geo_num,), dtype = np.float32)

    for i in range(ngeo):
        geo_x[i] = x
        geo_y[i] = y
        geo_z[i] = bottomz + i * voffset
        # print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' %(i, geo_x[i], geo_y[i], geo_z[i]))

    return geo_x, geo_y, geo_z

def plotarray(geo_x, geo_y, tit = 'Grid Array', xlen = 6000, ylen = 6000, offset = .0 ):

    from matplotlib import pyplot

    fig, ax = pyplot.subplots()
    ax.scatter(geo_x + offset, geo_y + offset, c='white', marker='v')
    ax.axis([0, xlen + 2*offset, 0, ylen + 2*offset])
    ax.patch.set_facecolor('#000080')
    pyplot.title(tit)
    pyplot.grid(True)
    pyplot.show()

def writearray(geo_x, geo_y, geo_z, geonum = 961, fname = "./gridarray.dat"):

    import codecs

    f = codecs.open(fname,'w','utf-8')
    for i in range(geonum):
        f.write(str(i) + '            ')
        f.write(str(geo_x[i])+'            ')
        f.write(str(geo_y[i]) + '            ')
        f.write(str(geo_z[i]) + '\r\n')
    f.close()

def wiggle(Data, tit, nt, ns, skipt=1, scale=0.5, ucolor='blue', lwidth=.1 ):
    """
	wiggle(Data,SH)
	"""

    import pylab
    from numpy import abs
    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    #
    max_val = np.zeros((nt,), dtype=np.float32)
    attn = np.zeros((nt,), dtype=np.float32)

    # calculate max for each trace
    for i in range(0, nt, skipt):
        trace = Data[:, i]
        if trace.max() > abs(trace.min()):
            max_val[i] = trace.max()
        else:
            max_val[i] = abs(trace.min())

    for i in range(0, nt, skipt):
        attn[i] = (max_val[i]-max_val.min())/(max_val.max()-max_val.min());

    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0

        coef = attn[i];

        trace = Data[:, i]
        if trace.max() > abs(trace.min()):
            maxval = trace.max()
        else:
            maxval = abs(trace.min())

        trace[0] = 0;
        trace[ns - 1] = 0
        pylab.plot(1 + i + trace / maxval * coef * scale, t, color=ucolor, linewidth=lwidth)
        for a in range(len(trace)):
            if (trace[a] < 0):
                trace[a] = 0;
            # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    pylab.title(tit)
    pylab.axis([0, nt + 1, ns, 1])

    pylab.grid(True)
    pylab.show()


def wiggleNorm(Data, tit, nt, ns, skipt=1, scale=0.5, ucolor='blue', lwidth=.1):
    """
	wiggle(Data,SH)
	"""

    import pylab
    from numpy import abs
    t = range(ns)
    # 	t = range(SH['ns'])*SH['dt']/1000000;
    for i in range(0, nt, skipt):
        #		trace=zeros(SH['ns']+2)
        #		dtrace=Data[:,i]
        #		trace[1:SH['ns']]=Data[:,i]
        #		trace[SH['ns']+1]=0
        trace = Data[:, i]
        if trace.max() > abs(trace.min()):
            maxval = trace.max()
        else:
            maxval = abs(trace.min())

        trace[0] = 0;
        trace[ns - 1] = 0
        pylab.plot(1 + i + trace / maxval * scale, t, color=ucolor, linewidth=lwidth)
        for a in range(len(trace)):
            if (trace[a] < 0):
                trace[a] = 0;
            # pylab.fill(i+Data[:,i]/maxval, t, 'k', linewidth=0)

    pylab.title(tit)
    pylab.axis([0, nt + 1, ns, 1])
    pylab.xlabel('Traces')
    pylab.ylabel('Samples')

    pylab.grid(True)
    # pylab.show()


def wavedisplay(tit = 'Waveform', fname = ' ', ntr = 961, ns = 3000, nt = 31, mode = 'norm'):

    import numpy as np
    from numpy import arange, reshape, transpose

    data = np.fromfile(fname, dtype=np.float32, count=-1, sep='')
    data_size = len(data)

    xcor = arange(1, ns + 1, 1)
    ylen = len(xcor)
    # print ylen

    seis = data.reshape(ntr, ns)
    vseis = transpose(seis)

    if mode == 'norm':
        wiggleNorm(vseis[:, 0:nt], tit, nt, ns, 1, 0.8, 'black', 0.5)

    else:
        wiggle(vseis[:, 0:nt], tit, nt, ns, 1, 0.8, 'black', 0.5)

## Read binary file containing seismic trace data


# wavedisplay('', '../output_data/gridarray_fwd3_syndata.dat', 961, 4000, 31)
# wavedisplay('', '../output_data/gridarray_fwd4_syndata.dat', 49, 4000, 7)

#Test Codes
# geox, geoy, geoz = dharray(12, 0, 0, 3000, 25)
# print(geox)


# ## paper code
# geo_x, geo_y, geo_z = gridarray(961, 6000, 6000)
# geox, geoy, geoz = gridarray(49, 6000, 6000)

# from matplotlib import pyplot
# offset = 10
# xlen = ylen = 6000
# fig1 = pyplot.subplot(221)
# fig1.scatter(geo_x + offset, geo_y + offset, c='white', marker='v')
# fig1.axis([0, xlen + 2*offset, 0, ylen + 2*offset])
# ax1 = pyplot.gca()
# ax1.set_xlabel('X(m)')
# ax1.set_ylabel('Y(m)')
# # ax.patch.set_facecolor('#000080')
# # pyplot.title(tit)
# # pyplot.grid(True)
# offset = 20
# fig2 = pyplot.subplot(223)
# fig2.scatter(geox + offset, geoy + offset, c='white', marker='v')
# fig2.axis([0, xlen + 2*offset, 0, ylen + 2*offset])
# ax2 = pyplot.gca()
# ax2.set_xlabel('X(m)')
# ax2.set_ylabel('Y(m)')

# fig3 = pyplot.subplot(222)
# wavedisplay('', '../output_data/gridarray_fwd3_syndata.dat', 961, 4000, 31)
# fig4 = pyplot.subplot(224)
# wavedisplay('', '../output_data/gridarray_fwd4_syndata.dat', 49, 4000, 7)

# pyplot.show()

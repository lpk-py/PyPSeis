'''

# A python script for calculating the probability function values
#  of microseismic location uncertainty
# 
# 
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

#Test codes
import numpy as np
from numpy import sqrt, square
from psmodules.pspicker import merpicker
from psmodules.pspdf import pdftp
from obspy.segy.segy import readSEGY
from psmodules.psarray import dharray

geonum = 12
gx, gy, gz = dharray(geonum, 0, 0, 3000, 25)


filename = '/Users/Uqer/Dropbox/SEGY/dr01.sgy'
st = readSEGY(filename)
tr=st.traces[0].data  
# pick = merpicker(tr, 512, 0.25, 20, 600, "true" )   

x1 = 500
x2 = 1000
z1 = 1000
z2 = 3000
vel = 3000

ns = 512
srate = 0.25
win = 20
threshold = 600
N = 12
to = 0
pSigma = 1
pdf = np.zeros((10000,), dtype = np.float32)
tp = np.zeros((geonum,), dtype = np.int32)
tmp = np.zeros((geonum,), dtype = np.int32)
ii = 0
jj = 0
for i in range(x1, x2, 5):
    for j in range(z1, z2, 20):
        for k in range(geonum):
            tp[k] = sqrt(square(i)+square(gz[k]-j))/vel*1000/srate
            # print('%d'%tp[k]) 
            tmp[k] = merpicker(st.traces[k*3].data, ns, srate, win, threshold, "false")
            # print('%d'%tmp[k]) 
        n = ii * 100 + jj
        print('%d'%n) 
        pdf[n] = pdftp(N, geonum, tp, tmp, to, pSigma)
        print('%100.99f'%pdf[n])  
        jj = jj + 1
    ii = ii + 1      
print(pdf)
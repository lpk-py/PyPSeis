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
from matplotlib import pyplot
import codecs



# Global parameters
X_LENGTH = 6000
Y_LENGTH = 6000
GEO_NUM = 1000
X_INNER = 3000
Y_INNER = 3000

## Generating Star Array
ARM_NUM = 12
RADIAN = 2*math.pi/ARM_NUM
OFFSET = 30
print ('Radian between arms is: %f' %(RADIAN))

RADIUS_INNER = 300
RADIUS_OUTER = 3000
NUM_GEO_PER_ARM = GEO_NUM/ARM_NUM
total_geo_num = NUM_GEO_PER_ARM * ARM_NUM
print('Number of geophones on each arm is: %d' %(NUM_GEO_PER_ARM))
print ('Total number of geophones used in Star-shaped Array is: %d' %(total_geo_num))

geo_x = np.zeros((total_geo_num,), dtype = np.float32)
geo_y = np.zeros((total_geo_num,), dtype = np.float32)
geo_z = np.zeros((total_geo_num,), dtype = np.float32)
for i in range(ARM_NUM):
    for j in range(NUM_GEO_PER_ARM):
        r = RADIUS_INNER + j * OFFSET
        geo_x[i * NUM_GEO_PER_ARM + j] = r * math.cos(i*RADIAN) + X_INNER
        geo_y[i * NUM_GEO_PER_ARM + j] = r * math.sin(i*RADIAN) + Y_INNER
        geo_z[i * NUM_GEO_PER_ARM + j] = 2.0
        # print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' % (i * NUM_GEO_PER_ARM + j, geo_x[i * NUM_GEO_PER_ARM + j] , geo_y[i * NUM_GEO_PER_ARM + j] , geo_z[i * NUM_GEO_PER_ARM + j] ))


fig, ax = pyplot.subplots()
ax.scatter(geo_x, geo_y, c='white', marker='v')
ax.axis([0, X_LENGTH, 0, Y_LENGTH])
ax.patch.set_facecolor('#0000AA')

print (geo_x)

f = codecs.open("../output_data/stararray.dat",'w','utf-8')
for i in range(total_geo_num):
    f.write(str(i) + '            ')
    f.write(str(geo_x[i])+'            ')
    f.write(str(geo_y[i]) + '            ')
    f.write(str(geo_z[i]) + '\r\n')
f.close()


## Generating Grid Array
# Case 1: Square Grid
ngeo = int(math.sqrt(GEO_NUM))
GEO_X_NUM = ngeo+1
GEO_Y_NUM = ngeo+1
print('Number of geophones in XY directions are [%d,%d]' %(GEO_X_NUM, GEO_Y_NUM))
X_OFFSET = X_LENGTH/GEO_X_NUM
print ('Offset in X direction is %d' %(X_OFFSET))
Y_OFFSET = Y_LENGTH/GEO_Y_NUM
print ('Offset in Y direction is %d' %(Y_OFFSET))
GLOBAL_OFFSET = 100.0
print('Geophone offsets in XY directions are [%d, %d]' %(X_OFFSET, Y_OFFSET))
total_geo_num = GEO_X_NUM*GEO_Y_NUM

geo_x = np.zeros((total_geo_num,), dtype = np.float32)
geo_y = np.zeros((total_geo_num,), dtype = np.float32)
geo_z = np.zeros((total_geo_num,), dtype = np.float32)

for i in range(GEO_Y_NUM):
    for j in range(GEO_X_NUM):
        if i>14 and i<18 and j>14 and j<18:
            geo_x[i * GEO_X_NUM + j] = 0.0
            geo_y[i * GEO_X_NUM + j] = 0.0
            geo_z[i * GEO_X_NUM + j] = 0.0
        else:
            geo_x[i * GEO_X_NUM + j] = j * X_OFFSET + GLOBAL_OFFSET
            geo_y[i * GEO_X_NUM + j] = i * Y_OFFSET + GLOBAL_OFFSET
            geo_z[i * GEO_X_NUM + j] = 2.0
        # print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' %(i * GEO_X_NUM + j, geo_x[i * GEO_X_NUM + j], geo_y[i * GEO_X_NUM + j], geo_z[i * GEO_X_NUM + j]))

NUM_GEO_REMOVED = 9
total_geo_num = GEO_X_NUM*GEO_Y_NUM - NUM_GEO_REMOVED
print ('Total number of geophones used in Square Array is: %d' %(total_geo_num))
geox = []
geoy = []
geoz = []
for i in range(GEO_X_NUM * GEO_Y_NUM):
    if geo_x[i] != 0:
        geox.append(geo_x[i])
    if geo_y[i] != 0:
        geoy.append(geo_y[i])
    if geo_z[i] != 0:
        geoz.append(geo_z[i])

# for i in range(len(geox)):
#     print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' % (i , geox[i], geoy[i], geoz[i]))


fig, ax = pyplot.subplots()
ax.scatter(geox, geoy, c='white', marker='v')
ax.axis([0, X_LENGTH, 0, Y_LENGTH])
ax.patch.set_facecolor('#000080')




# Case 2: Circular Grid

## Generating Patch Array
X_OFFSET = 25
Y_OFFSET = 22
GEO_X_NUM = 6
GEO_Y_NUM = 8
NUM_GEO_PER_PATCH = GEO_X_NUM * GEO_Y_NUM
RADIUS1 = 1000
NUM_PATCH_RADIUS1 = 4
RADIUS2 = 1800
NUM_PATCH_RADIUS2 = 7
RADIUS3 = 2700
NUM_PATCH_RADIUS3 = 10
PATCH_NUM = NUM_PATCH_RADIUS1 + NUM_PATCH_RADIUS2 + NUM_PATCH_RADIUS3
total_geo_num = NUM_GEO_PER_PATCH * PATCH_NUM
print ('Total number of geophones used in Patch Array is: %d' %(total_geo_num))

geox = []
geoy = []
geoz = []

RADIAN = 2*math.pi/NUM_PATCH_RADIUS1
for i in range(NUM_PATCH_RADIUS1):
    r = RADIUS1
    x0 = r * math.cos(i * RADIAN) + X_INNER
    y0 = r * math.sin(i * RADIAN) + Y_INNER
    for j in range(GEO_Y_NUM):
        for k in range(GEO_X_NUM):
            x = x0 + k * X_OFFSET
            geox.append(x)
            y = y0 + j * Y_OFFSET
            geoy.append(y)
            z = 2.0
            geoz.append(z)

RADIAN = 2*math.pi/NUM_PATCH_RADIUS2
for i in range(NUM_PATCH_RADIUS2):
    r = RADIUS2
    x0 = r * math.cos(i * RADIAN) + X_INNER
    y0 = r * math.sin(i * RADIAN) + Y_INNER
    for j in range(GEO_Y_NUM):
        for k in range(GEO_X_NUM):
            x = x0 + k * X_OFFSET
            geox.append(x)
            y = y0 + j * Y_OFFSET
            geoy.append(y)
            z = 2.0
            geoz.append(z)

RADIAN = 2*math.pi/NUM_PATCH_RADIUS3
for i in range(NUM_PATCH_RADIUS3):
    r = RADIUS3
    x0 = r * math.cos(i * RADIAN) + X_INNER
    y0 = r * math.sin(i * RADIAN) + Y_INNER
    for j in range(GEO_Y_NUM):
        for k in range(GEO_X_NUM):
            x = x0 + k * X_OFFSET
            geox.append(x)
            y = y0 + j * Y_OFFSET
            geoy.append(y)
            z = 2.0
            geoz.append(z)

# for i in range(len(geox)):
#     print('Geophone[%d] XYZ Coordinate: [%6.2f, %6.2f, %6.2f]' % (i , geox[i], geoy[i], geoz[i]))


fig, ax = pyplot.subplots()
ax.scatter(geox, geoy, c='white', marker='v')
ax.axis([0, X_LENGTH, 0, Y_LENGTH])
ax.patch.set_facecolor('#000080')

pyplot.show()

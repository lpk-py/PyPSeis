# !/Users/Uqer/anaconda/bin/python

# Test codes
import numpy as np
from math import floor, exp
from numpy import array, linspace, diff, append, \
    squeeze, concatenate, repeat, zeros, where, flipud, finfo
from psmodules.psraytrace import shooting
import matplotlib.pyplot as plt


#  Input geometry
xmin = 0
xmax = 2000
zmin = 0
zmax = 2000
# Make strata layer
zlayer = array([0, 80, 180, 300, 450, 630, 800,
                1050, 1200, 1400, 1620, 1820, 2000])
zlayer.shape = (len(zlayer), 1)
nlayer = len(zlayer)
layer = linspace(1, nlayer, nlayer)
# print(zlayer)
thick = abs(diff(zlayer))
x = append(xmin, xmax)
# print(x)
z = concatenate((zlayer, zlayer), axis=1)
# print(z)
dg = 10
ndg = (xmax - xmin) / dg + 1
# print(ndg)
xx = linspace(xmin, xmax, ndg)
nx = len(xx)
zz = repeat(zlayer, nx, axis=1)
# print(zz.shape)

# Source-Receiver Groups
# % Receiver Interval
dr = 20
# Front Spread Configuration
# Source
xs = array([0])
zs = array([0])
# # dxs = floor((xmax-xmin)/2)
# # ndxs = (xmax-xmin)/dxs + 1
# # xs = linspace(xmin, xmax, ndxs)
# # print(xs)
# zs = array([0, 0, 0])
ns = len(xs)
# Receiver
xr = linspace(100, 2000, 10)
# print(xr)
zr = zeros(len(xr))
nr = len(xr)
nray = ns * nr

# P wave velocity
vp = array([1800, 2000, 2200, 2500, 2400, 2700,
            3150, 2950, 3440, 3750, 4000, 4350, 600])
vp.shape = (len(vp), 1)
#  S wave velocity based on Castagna's rule
vs = (vp - 1360) / 1.16


#  Draw sources & receivers
fig, ax = plt.subplots()
cax = ax.pcolor(xx, zz, repeat(vp, nx, axis=1))
ax.plot(xs, zs, marker="*", color="black")
plt.hold(True)
ax.scatter(xr, zr, marker="v", color="red")
plt.hold(True)


# Run Ray Tracing

# Loop over for number of souce
for i in range(ns):
    # Loop over for number of receiver
    for j in range(nr):
        # % Loop over for number of layer
        for k in range(nlayer):
            print("k= %d" % k)
            #  Declare reflection boundary
            if zr[j] < zlayer[k] and zs[i] < zlayer[k]:
                zm = zz[k, :]
                # zf = min(zm) - finfo('float32').eps
                zf = min(zm) - exp(-10)
                # print("zf= %f" %zf)
                # Downgoing path
                ind = where(zlayer > zs[i])
                d = array(ind[0])
                d.shape = (len(d), 1)
                # print(d)
                if(len(d) == 0):
                    sdown = len(zlayer)
                else:
                    sdown = d[0]
                # print(sdown)
                ind = where(zlayer > zf)
                d = array(ind[0])
                d.shape = (len(d), 1)
                # print(d)
                if len(d) == 0:
                    edown = len(zlayer)
                else:
                    edown = d[0]
                # print(edown)

                if sdown + 1 > edown:
                    zd = zs[i]
                elif sdown + 1 == edown:
                    zd = append(zs[i], zlayer[sdown])
                    zd.shape = (len(zlayer[sdown]) + 1, 1)
                else:
                    zd = append(zs[i], zlayer[sdown:edown])
                    zd.shape = (len(zlayer[sdown:edown]) + 1, 1)
                # print(zd)
                nd = zd.size
                # print(nd)

                #  Upgoing path
                ind = where(zlayer > zr[j])
                u = array(ind[0])
                u.shape = (len(u), 1)
                # print(u)
                if(len(u) == 0):
                    sup = len(zlayer)
                else:
                    sup = u[0]
                # print(sup)
                ind = where(zlayer > zf)
                u = array(ind[0])
                u.shape = (len(u), 1)
                # print(u)
                if len(u) == 0:
                    eup = len(zlayer)
                else:
                    eup = u[0]
                # print(eup)

                if sup + 1 > eup + 1:
                    zu = zr[j]
                elif sup + 1 == eup + 1:
                    zu = append(zr[j], zlayer[sup])
                    zu.shape = (len(zlayer[sup]) + 1, 1)
                else:
                    zu = append(zr[j], zlayer[sup:eup + 1])
                    zu.shape = (len(zlayer[sup:eup + 1]) + 1, 1)
                # print(zu)
                nu = zu.size
                # print(nu)

                zn = array(append(zd, flipud(zu)))
                zn.shape = (len(zn), 1)
                # print(zn)

                #  Declare elastic parameter
                #  Downgoing elastic parameter
                if sdown == edown:
                    vpd = append(vp[sdown - 1], vp[edown - 1])

                else:
                    vpd = append(vp[sdown - 1:edown], vp[edown - 1])
                vpd.shape = (len(vpd), 1)
                # print(vpd)

                #  Upgoing elastic parameter
                if sup == eup:
                    vpu = array(append(vp[sup - 1], vp[eup - 1]))

                else:
                    vpu = array(append(vp[sup - 1:eup], vp[eup - 1]))
                vpu.shape = (len(vpu), 1)
                # print(vpu)

                #  Combine model elastic parameter
                vpp = array(append(vpd[0:len(vpd) - 1],
                                   flipud(vpu[0:len(vpu) - 1])))
                vpp.shape = (len(vpp), 1)
                # print(vpp)
                vps = vpp

                #  Start Raytracing (P-P, S-S, or P-S mode)
                ops = 1
                # ops=1 for PP mode; ops=2 for PS mode
                sx = array([xs[i]])
                # print(sx)
                rx = array([xr[j]])

                xh, zh, vh, pp, teta, time = shooting(
                    vpp, vps, zn, xx, sx, rx, ops)
                print(time)
                # theta = abs(teta); twt(k,j,i) = time;

                # Plot Ray
                if ops == 1:
                    plt.title("Seismic Raytracing (P-P mode)")
                    plt.plot(xh, zh)
                    plt.hold(True)


plt.axis([0, xmax, zmin - 0.03 * zmax, zmax])
plt.gca().invert_yaxis()
cbar = fig.colorbar(cax, orientation="horizontal")
cbar.set_label("Velocity (m/s)")

# plt.savefig("ref_seis_raytrace_model1.eps", format = "eps")
# plt.savefig("ref_seis_raytrace_model1.png", format = "png")
plt.show()


# !/Users/Uqer/anaconda/bin/python

# Test codes
import numpy as np
from math import floor, exp
from numpy import array, linspace, diff, append, \
    squeeze, concatenate, repeat, zeros, where, flipud, finfo
from psmodules.psraytrace import shooting, directshooting
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
xs = array([0, 0])
zs = array([1000, 2000])
# # dxs = floor((xmax-xmin)/2)
# # ndxs = (xmax-xmin)/dxs + 1
# # xs = linspace(xmin, xmax, ndxs)
# # print(xs)
# zs = array([0, 0, 0])
ns = len(xs)
# Receiver
xr = array([500, 1000])
# print(xr)
zr = array([0, 0])
nr = len(xr)
nray = ns * nr

# P wave velocity
vp = array([1800, 2000, 2200, 2500, 2400, 2700,
            3150, 2950, 3440, 3750, 4000, 4350, 4600])
vp.shape = (len(vp), 1)
#  S wave velocity based on Castagna's rule
vs = (vp - 1360) / 1.16


#  Draw sources & receivers
fig, ax = plt.subplots()
cax = ax.pcolor(xx, zz, repeat(vp, nx, axis=1))
ax.plot(xs, zs, markersize=10, marker="*", color="black")
plt.hold(True)
ax.scatter(xr, zr, s=50,  marker="v", color="red")
plt.axis([0, xmax, zmin - 0.03 * zmax, zmax])
plt.gca().invert_yaxis()
plt.hold(True)


# Run Ray Tracing

# Loop over for number of souce
for i in range(ns):
    # Loop over for number of receiver
    for j in range(nr):
        # Compare zs and zr to determine downgoing or upgoing shooting
        if zs[i] > zr[j]:
            #  Upgoing path
            ind = where(zlayer < zs[i])
            u = array(ind[0])
            u.shape = (len(u), 1)

            if len(u) == 0:
                eup = len(zlayer)
            else:
                eup = u[len(u) - 1]

            ind = where(zlayer > zr[j])
            u = array(ind[0])
            u.shape = (len(u), 1)

            if len(u) == 0:
                sup = len(zlayer)
            else:
                sup = u[0]

            if sup > eup:
                zu = zr[j]
            elif sup == eup:
                zuu = append(zr[j], zlayer[sup])
                zu = append(zu, zs[i])
                zu.shape = (3, 1)
            else:
                zuu = append(zr[j], zlayer[sup:eup])
                zu = append(zuu, zs[i])
                zu.shape = (len(zlayer[sup:eup]) + 2, 1)
            nu = len(zu)
            zn = flipud(zu)

            # Upgoing elastic parameter
            if sup - 1 == eup:
                vpu = vp[sup - 1]
                vsu = vs[sup - 1]
            else:
                vpu = vp[sup - 1:eup]
                vsu = vs[sup - 1:eup]

            # Combine model elastic parameter
            vpp = flipud(vpu)
            vps = flipud(vsu)

            # Start Raytracing (P-P, S-S, or P-S mode)
            ops = 1
            sx = array([xs[i]])
            rx = array([xr[j]])
            # ops=1 for PP mode; ops=2 for PS mode
            xh, zh, vh, pp, teta, time = shooting(
                vpp, vps, zn, xx, sx, rx,  ops)
            print(time)

        elif zs[i] == zr[j]:

            # Horizontal path
            ind = where(zlayer < zs[i])
            h = array(ind[0])
            h.shape = (len(h), 1)

            if len(h) == 0:
                hor = 1
            else:
                hor = h[len(h) - 1]

            zhor = append(zs[i], zr[j])
            nu = len(zhor)
            zn = zhor

            # Upgoing elastic parameter
            vph = vp[hor]
            vsh = vs[hor]

            # Combine model elastic parameter
            vpp = vph
            vps = vsh

            # Start Raytracing(P - P, S - S, or P - S mode)
            ops = 1
            # ops = 1 for PP mode, ops = 2 for PS mode
            sx = array([xs[i]])
            rx = array([xr[j]])
            xh, zh, vh, pp, teta, time = directshooting(
                vpp, vps, zn, xx, sx, rx, ops)
            print(time)

        else:
            # Downgoing path
            ind = where(zlayer > zs[i])
            d = array(ind[0])
            d.shape = (len(d), 1)

            if len(d) == 0:
                sdown = len(zlayer)
            else:
                sdown = d[0]

            ind = where(zlayer < zr[j])
            d = array(ind[0])
            d.shape = (len(d), 1)

            if len(d) == 0:
                edown = len(zlayer)
            else:
                edown = d[len(d) - 1]

            if sdown > edown:
                zd = zs[i]
            elif sdown == edown:
                zdd = append(zs[i], zlayer[sdown])
                zd = append(zdd, zr[j])
                zd.shape = (3, 1)
            else:
                zdd = append(zs[i], zlayer[sdown:edown])
                zd = append(zdd, zr[j])
                zd.shape = (len(zlayer[sdown:edown]) + 2, 1)
            nd = len(zd)
            zn = zd

            # Downgoing elastic parameter
            if sdown - 1 == edown:
                vpd = vp[sdown - 1]
                vsd = vs[sdown - 1]
            else:
                vpd = vp[sdown - 1: edown]
                vsd = vs[sdown - 1: edown]

            # Combine model elastic parameter
            vpp = vpd
            vps = vsd

            # Start Raytracing(P - P, S - S, or P - S mode)
            ops = 1
            # ops = 1 for PP mode, ops = 2 for PS mode
            sx = array([xs[i]])
            rx = array([xr[j]])
            xh, zh, vh, pp, teta, time = shooting(
                vpp, vps, zn, xx, sx, rx, ops)
            print(time)

        # Plot Ray
        if ops == 1:
            ax.set_title("Seismic Raytracing (P-P mode)")
            ax.plot(xh, zh, color="black")
            plt.hold(True)


cbar = fig.colorbar(cax, orientation="horizontal")
cbar.set_label("Velocity (m/s)")

# plt.savefig("ref_seis_raytrace_model1.eps", format = "eps")
# plt.savefig("ref_seis_raytrace_model1.png", format = "png")
plt.show()

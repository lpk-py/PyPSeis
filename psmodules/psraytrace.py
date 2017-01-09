'''

# A python module for seismic ray tracing
#
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
from math import exp, floor, ceil
from numpy.lib.scimath import sqrt
from numpy import array, linspace, ones, zeros, empty, repeat, \
    transpose, diff, where, real, cumsum, append, multiply, \
    arcsin, finfo, concatenate, square, flipud


def main():
    """
    # Executing main() function will run three examples, i.e.
    # 2D passive seismic rayracing example
    # 3D passive seismic raytracing example
    # 2D reflection Seismic raytracing example

    """

    # 3D passive seismic raytracing
    print("3D passive seismic raytracing example is running[Waiting...]")
    zlayer = array([0, 80, 180, 300, 450, 630, 800,
                    1050, 1200, 1400, 1620, 1820, 2000])
    dg = 10
    sourcex = array([3000, 3000])
    sourcey = array([3000, 3000])
    sourcez = array([1000, 2000])

    receiverx = array([300, 600])
    receivery = array([400, 800])
    receiverz = array([0, 0])

    # P wave velocity
    vp = array([1800, 2000, 2200, 2500, 2400, 2700,
                3150, 2950, 3440, 3750, 4000, 4350, 4600])

    #  S wave velocity based on Castagna's rule
    vs = (vp - 1360) / 1.16

    times = raytracing(vp, vs, zlayer, dg, sourcex, sourcey,
                       sourcez, receiverx, receivery, receiverz)
    print("One Way Travel Time is :")
    print(times)
    print("3D passive seismic raytracing example running completed[OK]")

    # 2D passive seismic raytracing
    print("2D passive seismic raytracing example is running[Waiting...]")

    # Make strata layer
    zlayer = array([0, 80, 180, 300, 450, 630, 800,
                    1050, 1200, 1400, 1620, 1820, 2000])
    zlayer.shape = (len(zlayer), 1)
    nlayer = len(zlayer)
    layer = linspace(1, nlayer, nlayer)

    #  Input geometry
    xmin = 0
    xmax = 6000
    zmin = 0
    zmax = 2000

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
    xs = array([3000])
    zs = array([1000])
    # # dxs = floor((xmax-xmin)/2)
    # # ndxs = (xmax-xmin)/dxs + 1
    # # xs = linspace(xmin, xmax, ndxs)
    # # print(xs)
    # zs = array([0, 0, 0])
    ns = len(xs)
    # Receiver
    xr = array([500])
    # print(xr)
    zr = array([0])
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
                    zuu = append(zr[j], zlayer[sup:eup + 1])
                    zu = append(zuu, zs[i])
                    zu.shape = (len(zlayer[sup:eup + 1]) + 2, 1)
                nu = len(zu)
                zn = flipud(zu)

                # Upgoing elastic parameter
                if sup - 1 == eup:
                    vpu = vp[sup - 1]
                    vsu = vs[sup - 1]
                else:
                    vpu = vp[sup - 1:eup + 1]
                    vsu = vs[sup - 1:eup + 1]

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
    print("2D passive seismic raytracing example running completed[OK]")

    # 2D reflection seismic raytracing
    print("2D reflection seismic raytracing example is running[Waiting...]")
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
                # print("k= %d" % k)
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
    print("2D reflection seismic raytracing example running completed[OK]")


# Functions
def raytracing(vp, vs, zlayer, dg, sourcex, sourcey, sourcez, receiverx, receivery, receiverz):

    #  Input geometry
    xmin = 0
    xmax = 20000
    zmin = 0
    zmax = zlayer.max()
    ndg = (xmax - xmin) / dg + 1

    # Make strata layer
    zlayer.shape = (len(zlayer), 1)
    nlayer = len(zlayer)
    layer = linspace(1, nlayer, nlayer)

    thick = abs(diff(zlayer))
    x = append(xmin, xmax)
    z = concatenate((zlayer, zlayer), axis=1)

    xx = linspace(xmin, xmax, ndg)
    nx = len(xx)
    zz = repeat(zlayer, nx, axis=1)

    # P wave velocity
    vp.shape = (len(vp), 1)
    #  S wave velocity based on Castagna's rule
    # vs = (vp - 1360) / 1.16
    vs.shape = (len(vs), 1)

    # Source-Receiver Groups
    # Source
    zs = sourcez
    xs = sourcex - sourcex
    ns = len(xs)
    # Receiver
    zr = receiverz
    # print(zr)
    nr = len(zr)
    xr = empty((nr, 1), dtype="float32")

    nray = ns * nr

    times = []
    thetas = []

    # Run Ray tracing
    # # Loop over for number of souce
    for i in range(ns):
        # Loop over for number of receiver
        for j in range(nr):
            xr[j] = sqrt((sourcex[i] - receiverx[j]) * (sourcex[i] - receiverx[j]) +
                         (sourcey[i] - receivery[j]) * (sourcey[i] - receivery[j]))
            # print(xr[j])
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
                xh, zh, vh, pp, teta, ttime = shooting(
                    vpp, vps, zn, xx, sx, rx,  ops)
                if xs[i] == xr[j]:
                    zv = abs(diff(zn, axis=0))
                    tv = zv / vpp
                    tt = cumsum(tv)
                    tt_size = tt.size
                    ttime = tt[tt_size - 1]

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
                xh, zh, vh, pp, teta, ttime = directshooting(
                    vpp, vps, zn, xx, sx, rx, ops)

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
                xh, zh, vh, pp, teta, ttime = shooting(
                    vpp, vps, zn, xx, sx, rx, ops)
                if xs[i] == xr[j]:
                    zv = abs(diff(zn, axis=0))
                    tv = zv / vpp
                    tt = cumsum(tv)
                    tt_size = tt.size
                    ttime = tt[tt_size - 1]

            # Store traveltimes and incidence angles
            # print(time)
            times.append(ttime)
            t = array(times, dtype='float32')
            t.shape = (len(times), 1)
            thetas.append(abs(teta[len(teta) - 1]))
            teta = array(thetas, dtype='float32')
            teta.shape = (len(thetas), 1)

    return t, teta


def directshooting(vpp, vps, zn, xx, xs, xr, ops):
    # Horizontal path
    if xs < xr:
        xh = [xs, xr]
    else:
        xh = [xr, xs]

    zh = zn
    vh = vpp
    teta = 0.0

    if ops == 1:
        pp = 1 / vpp
        time = abs(xs - xr) / vpp
    else:
        pp = 1 / vps
        time = abs(xs - xr) / vps

    return xh, zh, vh, pp, teta, time


def shooting(vpp, vps, zn, xx, xs, xr, ops):
    # some constants
    itermax = 50
    offset = abs(xs - xr)
    # print("offset:")
    # print(len(offset))
    xc = 10

    # determin option
    if (ops == 1):
        vh = vpp
    elif (ops == 2):
        vh = vps

    # initial guess of the depth & time
    zh = zn - finfo("float32").eps
    # zh = zn - exp(-10)
    # print("zh:")
    # print(zh)
    t = float("inf") * ones((len(offset),), dtype=np.float32)
    # t = exp(100) * ones((len(offset),), dtype=np.float32)
    # print (t)
    p = float("inf") * ones((len(offset),), dtype=np.float32)
    # p = exp(100) * ones((len(offset),), dtype=np.float32)

    # start raytracing
    # trial shooting
    pmax = 1 / min(vh)
    # print("pmax= %10.9f" %(pmax))
    pp = np.linspace(0, 1 / max(vh), len(xx))
    pp.shape = (1, len(pp))
    # print("pp: ")
    # print (pp)

    sln = vh[0:len(zh)] * pp  # - exp(-20)
    # print("sln: ")
    # print(sln)
    vel = vh[0:len(zh) - 1] * ones((1, len(pp)))
    # print ("vel: ")
    # print(vel)
    dz = abs(diff(zh, axis=0)) * ones((1, len(pp)))
    # print("dz:")
    # print (dz.size)
    dim_sln = sln.shape
    # print("dim_sln:")
    # print (dim_sln)
    if (dim_sln[0] > 1):
        xn = sum((dz * sln) / sqrt(1 - sln**2))
        # print("xn:")
        # print(xn)
        tt = sum(dz / (vel * sqrt(1 - sln**2)))
        # print("tt:")
        # print(tt)
    else:
        xn = (dz * sln) / sqrt(1 - sln**2)
        tt = dz / (vel * sqrt(1 - sln**2))
    xn.shape = (1, len(xn))
    tt.shape = (1, len(tt))

    # print(xn.shape)
    xmax = xn.max()
    # print("xmax= %f" %(xmax))

    # bisection method
    # start bisection method

    for k in range(len(offset)):
        # print("k= %d" %k)
        # analyze the radius of target
        n = xn.size
        # print("n= %d" %n)
        xa = xn[:, 0:n - 1]
        # print("xa:")
        # print(xa.shape)
        xb = xn[:, 1:n]
        # print("xb:")
        # print(xb)
        opt1 = empty((1, n - 1))
        opt2 = empty((1, n - 1))
        opt = empty((1, n - 1))
        for i in range(n - 1):
            if xa[:, i] <= offset[k] and xb[:, i] > offset[k]:
                opt1[:, i] = 1
            else:
                opt1[:, i] = 0

            if xa[:, i] >= offset[k] and xb[:, i] < offset[k]:
                opt2[:, i] = 1
            else:
                opt2[:, i] = 0

        opt = opt1 + opt2
        # print("opt1:")
        # print(opt1)
        # print("opt2:")
        # print(opt2)
        # print("opt:")
        # print(opt.shape)

        index = where(opt == 1)
        ind = index[1]
        # print (ind)
        # print("ind= %d" %ind)
        if len(index) == 0:
            if (offset(k) >= xmax):
                a = n
                b = []
            else:
                a = []
                b = 1
        else:
            a = ind
            b = ind + 1
        # print(type(a))
        x1 = xn[:, a][0, 0]
        x2 = xn[:, b][0, 0]
        t1 = tt[:, a][0, 0]
        t2 = tt[:, b][0, 0]
        p1 = pp[:, a][0, 0]
        p2 = pp[:, b][0, 0]
        iter = 0
        err = (b - a) / 2
        # print("err= %f" %err)

        # while((iter < itermax) and (abs(err)<1)):
        while((iter < itermax) and abs(err) < 1):
            # print(iter)
            iter = iter + 1
            xt1 = abs(offset[k] - x1)
            xt2 = abs(offset[k] - x2)
            if (xt1 < xc) and (xt1 <= xt2):
                # linear interpolation
                t[k] = t1 + (offset[k] - x1) * (t2 - t1) / (x2 - x1)
                p[k] = p1 + (offset[k] - x1) * (p2 - p1) / (x2 - x1)
            elif (xt2 < xc) and (xt2 <= xt1):
                t[k] = t2 + (offset[k] - x2) * (t1 - t2) / (x1 - x2)
                p[k] = p2 + (offset[k] - x2) * (p1 - p2) / (x1 - x2)
            # set new ray parameter
            if a.size == 0:
                p2 = p1
                p1 = 0
            elif b.size == 0:
                p1 = p2
                p2 = pmax
            # print ("p1= %f  p2= %f" %(p1, p2))
            pnew = linspace(array([p1, p2]).min(), array([p1, p2]).max(), 3)
            pnew.shape = (1, len(pnew))
            # print("pnew:")
            # print (pnew)

            # do shooting by new ray parameter
            sln = vh[0:len(zh)] * pnew[:, 1]
            # print("sln:")
            # print(sln)
            vel = vh[0:len(zh)] * ones((1, 1))
            # print("vel:")
            # print(vel)
            dz = abs(diff(zh, axis=0)) * ones((1, 1))
            # print("dz:")
            # print (dz)
            dim_sln = sln.shape
            # print("dim_sln:")
            # print (dim_sln)
            if (dim_sln[0] > 1):
                xtemp = sum((dz * sln) / sqrt(1 - sln**2))
                # print("xn:")
                # print(xn)
                ttemp = sum(dz / (vel * sqrt(1 - sln**2)))
                # print("tt:")
                # print(tt)
            else:
                xtemp = (dz * sln) / sqrt(1 - sln**2)
                ttemp = dz / (vel * sqrt(1 - sln**2))
            xnew = array([x1, xtemp, x2])
            xnew.shape = (1, len(xnew))

            # print("xnew: ")
            # print(xnew)
            tnew = array([t1, ttemp, t2])
            tnew.shape = (1, len(tnew))
            # tnew = [t1, ttemp, t2]
            # print(tnew)
            xmax = xnew.max()
            # print("xmax= %f" %xmax)

            # analyze the radius of target
            n = xnew.size
            # print("n= %d" %n)
            xa = xnew[:, 0:n - 1]
            # print("xa:")
            # print(xa)
            xb = xnew[:, 1:n]
            # print("xb:")
            # print(xb)
            opt1 = empty((1, n - 1))
            opt2 = empty((1, n - 1))
            opt = empty((1, n - 1))

            for i in range(n - 1):
                # print("i= %d" %i)
                # print(offset[k])
                if xa[:, i] <= offset[k] and xb[:, i] > offset[k]:
                    opt1[:, i] = 1
                else:
                    opt1[:, i] = 0
                if xa[:, i] >= offset[k] and xb[:, i] < offset[k]:
                    opt2[:, i] = 1
                else:
                    opt2[:, i] = 0

            opt = opt1 + opt2
            # print("opt1:")
            # print(opt1)
            # print("opt2:")
            # print(opt2)
            # print("opt:")
            # print(opt.shape)

            index = where(opt == 1)
            ind = index[1]
            # print (ind)

            a = ind
            b = ind + 1
            x1 = xnew[:, a][0, 0]
            x2 = xnew[:, b][0, 0]
            # print(tnew)
            t1 = tnew[:, a][0, 0]
            t2 = tnew[:, b][0, 0]
            p1 = pnew[:, a][0, 0]
            p2 = pnew[:, b][0, 0]
            err = (b - a) / 2

            # declare ray parameter
            if xr[0] > xs[0]:
                pp = p
            else:
                pp = -p
            # compute travel time & angle
            dx = real((pp * vh * dz) / sqrt(1 - pp * pp * vh * vh))
            # print("dx:")
            # print(dx)
            xx = xs + cumsum(dx)
            # print("xx:")
            # print(xx)
            xh = (append(xs, xx)).reshape(xs.size + xx.size, 1)
            # print ("xh:")
            # print(xh)
            dz = real(dx * sqrt(1 - pp * pp * vh * vh)) / (pp * vh)
            dt = dz / (vh * sqrt(1 - pp * pp * vh * vh))
            tt = cumsum(dt)
            tt_size = tt.size
            time = tt[tt_size - 1]
            # print("time= %f" %time)

            teta = real(arcsin(multiply(pp, vh)))
            # print("teta:")
            # print(teta)

    return xh, zh, vh, pp, teta, time


# This will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Author: Thomas Meier

import numpy as np
import matplotlib.pyplot as plt
from scalingLawsReinhardt import *
import argparse

if __name__ == '__main__':
    parameters = initParameters()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-gamma', nargs='?', type=float, default = 1.0, const = 1.0, help='Impactor/target mass fraction, default=1.0')
    parser.add_argument('-Vimin', nargs='?', type=float, default = 20, const = 20, help='minimum impact velcity [km/s], default=20')
    parser.add_argument('-Vimax', nargs='?', type=float, default = 100, const = 100, help='maximum impact velocity [km/s], default=100')
    parser.add_argument('-Mtargmin', nargs='?', type=float, default = 0.1, const = 0.1, help='minimum target mass [ME], default=0.1')
    parser.add_argument('-Mtargmax', nargs='?', type=float, default = 10, const = 10, help='maximum target mass [ME], default=10')
    parser.add_argument('-NVi', nargs='?', type=int, default = 201, const = 201, help='number of points along the Vi axis, default=201')
    parser.add_argument('-NMtarg', nargs='?', type=int, default = 201, const = 201, help='number of points along the Mtarg axis, default=201')
    parser.add_argument('-plot', dest='plot', default=False, action='store_true', help='adds a plot, but needs more time to compute')
    args = parser.parse_args()

    # Apply arguments
    gamma = args.gamma
    Vimin = args.Vimin
    Vimax = args.Vimax
    Mtargmin = args.Mtargmin
    Mtargmax = args.Mtargmax
    NVi = args.NVi
    NMtarg = args.NMtarg
    plot = args.plot
    
    # Vectors
    Vis = np.linspace(Vimin,Vimax,NVi)
    Mtargs = np.logspace(np.log10(Mtargmin),np.log10(Mtargmax),NMtarg)
    
    # Results
    Zs = np.zeros(len(Mtargs)*len(Vis))
    Mlrgs = np.zeros(len(Mtargs)*len(Vis))
    flag1 = np.zeros(len(Mtargs)*len(Vis))
    flag2 = np.zeros(len(Mtargs)*len(Vis))
    flag3 = np.zeros(len(Mtargs)*len(Vis))
    flag4 = np.zeros(len(Mtargs)*len(Vis))
    qq = -1;
    for i in range(len(Mtargs)):
        for j in range(len(Vis)):
            qq = qq + 1
            if (i == 0) :
                flag1[qq] = 1
            if (j == 0) :
                flag2[qq] = 1
            if (i == len(Mtargs)-1) :
                flag3[qq] = 1
            if (j == len(Vis)-1) :
                flag4[qq] = 1
            if (plot or flag1[qq] or flag2[qq] or flag3[qq] or flag4[qq]):
                QRbyQRDstar = calcQRbyQRDstar(Mtargs[i],gamma,Vis[j],parameters)
                Mlrgs[qq] = calcMlrgofMtarggammaVi(Mtargs[i],gamma,Vis[j],parameters)
                Zs[qq] = ZofQRbyQRDstar(QRbyQRDstar,parameters)
            
    upperM = Mlrgs[((flag1 == 1) | (flag2 == 1))]
    upperZ = Zs[((flag1 == 1) | (flag2 == 1))]
    lowerM = Mlrgs[((flag3 == 1) | (flag4 == 1))]
    lowerZ = Zs[((flag3 == 1) | (flag4 == 1))]
    
    ind = np.argsort(upperM)
    upperM = upperM[ind]
    upperZ = upperZ[ind]
    ind = np.argsort(lowerM)
    lowerM = lowerM[ind]
    lowerZ = lowerZ[ind]

    f = open('lower.txt','w')
    f.write('Mlarg Z\n')
    for i in range(len(lowerM)):
        f.write('{:.15e} {:.15e}\n'.format(lowerM[i],lowerZ[i]))
    f.close()
    f = open('upper.txt','w')
    f.write('Mlarg Z\n')
    for i in range(len(upperM)):
        f.write('{:.15e} {:.15e}\n'.format(upperM[i],upperZ[i]))
    f.close()
    
    # plot
    if plot:
        plt.figure()
        plt.scatter(Mlrgs,Zs,1,label='Sampling points')
        plt.fill(np.append(lowerM, upperM[::-1]), np.append(lowerZ, upperZ[::-1]), color='grey', linewidth=0.0, alpha=0.6, zorder=10,label='Possible region')
        plt.legend()
        plt.xlabel(r"Fragment mass $\left[ \mathrm{M_{\oplus}} \right]$")
        plt.ylabel(r"Iron mass fraction $\left[ \mathrm{M_{Fe}/M_{frg}} \right]$")
        plt.gca().set_xlim(0.0,np.maximum(np.amax(lowerM),np.amax(upperM)))
        plt.gca().set_ylim(0.3,1.0)
        plt.show()
    
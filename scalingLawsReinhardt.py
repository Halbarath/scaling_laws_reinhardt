#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Author: Thomas Meier

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def calcQRDstar(Mtarg,gamma,Vi,parameters):
    RC1 = (3 * (1 + gamma) * Mtarg / (4 * np.pi * parameters['rho0'])) ** (1.0 / 3.0)
    QRDstar = parameters['qg'] * RC1 ** (3 * parameters['mu']) * Vi ** (2 - 3 * parameters['mu'])
    return QRDstar

def calcQR(gamma,Vi):
    QR = 1.0 / 2.0 * gamma / (1 + gamma) ** 2 * Vi **2
    return QR

def calcQRbyQRDstar(Mtarg,gamma,Vi,parameters):
    QR = calcQR(gamma,Vi)
    QRDstar = calcQRDstar(Mtarg,gamma,Vi,parameters)
    QRbyQRDstar = QR / QRDstar
    return QRbyQRDstar

def ZofQRbyQRDstar(QRbyQRDstar,parameters):
    aFesc = 1 - (parameters['ZFei'] + parameters['aFe'] * parameters['qsc'] ** parameters['bFe'])
    if (QRbyQRDstar < parameters['qsc']):
        Z = parameters['ZFei'] + parameters['aFe'] * QRbyQRDstar ** parameters['bFe']
    else:
        Z = 1 - aFesc * np.exp(-parameters['bFEsc'] * (QRbyQRDstar - parameters['qsc']))
    return Z

def calcMlrgofMtarggammaVi(Mtarg,gamma,Vi,parameters):
    QRbyQRDstar = calcQRbyQRDstar(Mtarg,gamma,Vi,parameters)
    if (QRbyQRDstar < parameters['b']):
        beta = 0.5 * np.cos(QRbyQRDstar * np.pi / 2.0) + 0.5;
    else:
        beta = (0.5 * np.cos(parameters['b'] * np.pi / 2.0) + 0.5) / (parameters['b'] ** parameters['c']) * QRbyQRDstar ** parameters['c']
    Mlrg = ((1+gamma)*beta)* Mtarg
    return Mlrg

def calcMtargofMlrggammaVi(Mtarg,Mlrg,gamma,Vi,parameters):
    QRbyQRDstar = calcQRbyQRDstar(Mtarg,gamma,Vi,parameters)
    if (QRbyQRDstar < parameters['b']):
        beta = 0.5 * np.cos(QRbyQRDstar * np.pi / 2.0) + 0.5;
    else:
        beta = (0.5 * np.cos(parameters['b'] * np.pi / 2.0) + 0.5) / (parameters['b'] ** parameters['c']) * QRbyQRDstar ** parameters['c']
    Mtarg = Mlrg / ((1 + gamma) * beta)
    return Mtarg

def ZofMlrggammaVi(Mlrg,gamma,Vi,parameters):
    f = lambda x: calcMtargofMlrggammaVi(x,Mlrg,gamma,Vi,parameters)-x
    sol = root_scalar(f,method='brentq',bracket=[1e-6,1e6])
    Mtarg = sol.root
    QRbyQRDstar = calcQRbyQRDstar(Mtarg,gamma,Vi,parameters)
    Z = ZofQRbyQRDstar(QRbyQRDstar,parameters)
    return Z, Mtarg

def ViofMlrggammaZ(Mlrg,gamma,Z,parameters):
    f = lambda x: ZofMlrggammaVi(Mlrg,gamma,x,parameters)[0]-Z
    sol = root_scalar(f,method='brentq',bracket=[0.01,1e3])
    Vi = sol.root
    return Vi

def initParameters():
    # Parameters
    parameters = {}
    # Parameters for Mtarg(Mlrg,gamma,Vi)
    parameters['b'] = 1.2645
    parameters['c'] = -8.1214
    # Parameters for QRD*(Mtarg,gamma,Vi)
    parameters['rho0'] = 5.972e-12
    parameters['qg'] = 1.3299e-05 # this is in our unit system, paper has it in cgs
    parameters['mu'] = 0.6164
    # Parameters for Z(QR/QRD*)
    parameters['ZFei'] = 0.3
    parameters['aFe'] = 0.2099
    parameters['bFe'] = 2.1516
    parameters['qsc'] = 1.2645
    parameters['bFEsc'] = 6.4890
    return parameters


if __name__ == '__main__':
    parameters = initParameters()

    # figure 4
    plt.figure(1)
    MlrgMin = 1e-4
    MlrgMax = 3e2
    Mlrgs = np.logspace(np.log10(MlrgMin),np.log10(MlrgMax),101)
    gamma = 1.0
    Zs = np.array([0.301,0.4,0.5,0.7,0.8,0.9,0.95,0.99])
    Vis = np.zeros((len(Zs),len(Mlrgs)))
    Mtarg = np.zeros((len(Zs),len(Mlrgs)))
    for j in range(len(Zs)):
        for i in range(len(Mlrgs)):
            Vis[j][i] = ViofMlrggammaZ(Mlrgs[i],gamma,Zs[j],parameters)
            Mtarg[j][i] = ZofMlrggammaVi(Mlrgs[i],gamma,Vis[j][i],parameters)[1]

    for j in range(len(Zs)):
            x = Mlrgs
            y = Vis[j][:]
            plt.plot(x,y,label='Z = {}'.format(Zs[j]))

    plt.xscale('log')
    plt.gca().set_xlim(MlrgMin,MlrgMax)
    plt.xlabel(r"Fragment mass $\left[ \mathrm{M_{\oplus}} \right]$")
    plt.ylabel(r"Impact velocity $\left[ \mathrm{km\ s^{-1}} \right]$")
    plt.legend()

    # figure 5 and 6
    plt.figure(2)
    plt.figure(3)
    MlrgMin = 1e-4
    MlrgMax = 3e2
    Mlrgs = np.logspace(np.log10(MlrgMin),np.log10(MlrgMax),101)
    gammas = np.array([0.1,0.25,1.0])
    Vis = np.array([20,40,60,80,160])
    Zs = np.zeros((len(Vis),len(gammas),len(Mlrgs)))
    Mtarg = np.zeros((len(Vis),len(gammas),len(Mlrgs)))
    for k in range(len(gammas)):
        for j in range(len(Vis)):
            for i in range(len(Mlrgs)):
                Zs[j][k][i],Mtarg[j][k][i] = ZofMlrggammaVi(Mlrgs[i],gammas[k],Vis[j],parameters)

    for k in range(len(gammas)):
        for j in range(len(Vis)):
            x = Mlrgs
            y = Zs[j][k][:]
            plt.figure(2)
            plt.plot(x,y,label='Vi = {}, gamma = {}'.format(Vis[j],gammas[k]))
            y = Mtarg[j][k][:]
            plt.figure(3)
            plt.plot(x,y,label='Vi = {}, gamma = {}'.format(Vis[j],gammas[k]))
    plt.figure(2)
    plt.xscale('log')
    plt.gca().set_xlim(MlrgMin,MlrgMax)
    plt.gca().set_ylim(0.3,1)
    plt.legend()
    plt.xlabel(r"Fragment mass $\left[ \mathrm{M_{\oplus}} \right]$")
    plt.ylabel(r"Iron mass fraction $\left[ \mathrm{M_{Fe}/M_{frg}} \right]$")
    plt.figure(3)
    plt.xscale('log')
    plt.gca().set_xlim(MlrgMin,MlrgMax)
    plt.legend()
    plt.xlabel(r"Fragment mass $\left[ \mathrm{M_{\oplus}} \right]$")
    plt.ylabel(r"Target mass $\left[ \mathrm{M_{\oplus}} \right]$")
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
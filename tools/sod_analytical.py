# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 12:03:03 2022

@author: Leo Basov
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

def _calc_u3(p1, p3, rhoL, Gamma, gamma, beta):
    num = (1 - Gamma**2) * p1**(1.0/gamma)
    denum = Gamma**2 * rhoL
    return (p1**beta - p3**beta) * np.sqrt(num / denum)

def _calc_u4(p3, p5, rhoR, Gamma):
    return (p3 - p5) * np.sqrt((1 - Gamma) / (rhoR*(p3 - Gamma*p5)))

def _calc_p3(p, rho, Gamma, gamma, beta):
    p1 = p[0]
    p3 = p[0]
    p5 = p[4]
    dp = 1e-5
    rhoL = rho[0]
    rhoR = rho[-1]
    err = 1e+5
    err_last = 1e+6
    
    while err < err_last:
        err_last = err
        err = abs(_calc_u3(p1, p3, rhoL, Gamma, gamma, beta) - _calc_u4(p3, p5, rhoR, Gamma))
        p3 -= dp
     
    return p3

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('t', type=float, help='time of results')

    args = parser.parse_args()
    
    # number points per segment
    N = 10
    
    # gas = Ar
    gamma = 1.67
    Gamma = (gamma - 1.0) / (gamma + 1.0)
    beta = (gamma - 1.0) / (2.0 * gamma)
    
    # boundary conditions
    rho = np.zeros(5)
    p = np.zeros(5)
    u = np.zeros(5)
    
    rho[0] = 1.0
    p[0] = 1.0
    u[0] = 0.0
    
    rho[-1] = 0.125
    p[-1] = 0.1
    u[-1] = 0.0
    
    # calculating states
    p[2] = _calc_p3(p, rho, Gamma, gamma, beta)
    p[3] = p[2]
    
    # speed of characterisics
    u[1] = np.sqrt(gamma * p[0] / rho[0]) # speed of rarefication going left
    u[2] = 0
    u[3] = (p[2] - p[4]) / np.sqrt((rho[4]/2) * ((gamma + 1) * p[2] + (gamma - 1) * p[4]))
    u[4] = np.sqrt(gamma * p[-1] / rho[-1]) # speed of pressure increase going right
    
    rho[2] = rho[0]*(p[2] / p[0])**(1.0/gamma)
    rho[3] = rho[4] * (p[3] + Gamma*p[4]) / (p[4] + Gamma*p[3])
    
    # calc x
    x = np.array([0.0, 0.5 - args.t*u[1], 0.5, 0.5 + args.t*u[3], 0.5 + args.t*u[4], 1,0])
    
    # calc rho
    rho1 = np.ones(N) * rho[0]
    rho2 = np.linspace(rho[0], rho[2], N)
    rho3 = np.ones(N) * rho[2]
    rho4 = np.ones(N) * rho[3]
    rho5 = np.ones(N) * rho[4]
    
    
    plt.plot(np.linspace(x[0], x[1], N), rho1)
    plt.plot(np.linspace(x[1], x[2], N), rho2)
    plt.plot(np.linspace(x[2], x[3], N), rho3)
    plt.plot(np.linspace(x[3], x[4], N), rho4)
    plt.plot(np.linspace(x[4], x[5], N), rho5)
    
    plt.show()
    
    print("done")
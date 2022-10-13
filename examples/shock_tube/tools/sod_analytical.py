# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 12:03:03 2022

@author: Leo Basov
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

kb = 1.380649e-23

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
    dp = p1 * 1e-3
    rhoL = rho[0]
    rhoR = rho[-1]
    err = 1e+5
    err_last = 1e+6
    
    while err < err_last:
        err_last = err
        err = abs(_calc_u3(p1, p3, rhoL, Gamma, gamma, beta) - _calc_u4(p3, p5, rhoR, Gamma))
        p3 -= dp
     
    return p3

def plot_val(x, val, name):
    plt.plot((x[0], x[1]), (val[0], val[1]))
    plt.plot((x[1], x[2]), (val[1], val[2]))
    plt.plot((x[2], x[3]), (val[2], val[2]))
    plt.plot((x[3], x[3]), (val[2], val[3]))
    plt.plot((x[3], x[4]), (val[3], val[3]))
    plt.plot((x[4], x[4]), (val[3], val[4]))
    plt.plot((x[4], x[5]), (val[4], val[4]))
    
    plt.ylabel(name)
    plt.xlabel("x")
    
    plt.show()
    
def write_values(x, val, name, file):
    with open(file + "_" + name + ".csv", "w") as file:
        file.write("{}, {}\n".format(x[0], val[0]))
        file.write("{}, {}\n".format(x[1], val[1]))
        file.write("{}, {}\n".format(x[2], val[2]))
        file.write("{}, {}\n".format(x[3], val[2]))
        file.write("{}, {}\n".format(x[3], val[3]))
        file.write("{}, {}\n".format(x[4], val[3]))
        file.write("{}, {}\n".format(x[4], val[4]))
        file.write("{}, {}\n".format(x[5], val[4]))
    
def check_args(args):
    if args.p is not None and (args.p[0] <= args.p[1]):
        raise Exception("p1 must be > p2")
    elif args.rho is not None and (args.rho[0] <= args.rho[1]):
        raise Exception("rho1 must be > rho2")
    elif args.L is not None and (args.L <= 0.0):
        raise Exception("L must be > 0")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=float, help='time of results', default=3.0e-5)
    parser.add_argument('-p', type=float, help='pressure', nargs = 2, default=(100, 10.0))
    parser.add_argument('-rho', type=float, help='density', nargs = 2, default=(0.0016036396304, 0.0016036396304*0.1))
    parser.add_argument('-L', type=float, help='tube length', default=0.1)
    parser.add_argument('-plt', type=bool, help='plot values', const=True, nargs='?')
    parser.add_argument('-w', type=str, help='write to file')

    args = parser.parse_args()
    
    check_args(args)
    
    # number points per segment
    N = 10
    
    # gas = Ar
    gamma = 1.67
    Gamma = (gamma - 1.0) / (gamma + 1.0)
    beta = (gamma - 1.0) / (2.0 * gamma)
    mass = 6.6422e-26
    
    # boundary conditions
    rho = np.zeros(5)
    p = np.zeros(5)
    u = np.zeros(5)
    
    rho[0] = args.rho[0]
    p[0] = args.p[0]
    u[0] = 0.0
    
    rho[-1] = args.rho[1]
    p[-1] = args.p[1]
    u[-1] = 0.0
    
    # calculating states
    p[1] = p[0]
    p[2] = _calc_p3(p, rho, Gamma, gamma, beta)
    p[3] = p[2]
    
    # speed of characterisics
    u[1] = np.sqrt(gamma * p[0] / rho[0]) # speed of rarefication going left
    u[2] = 0
    u[3] = (p[2] - p[4]) / np.sqrt((rho[4]/2) * ((gamma + 1) * p[2] + (gamma - 1) * p[4]))
    u[4] = np.sqrt(gamma * p[-1] / rho[-1]) # speed of pressure increase going right
    
    rho[1] = rho[0]
    rho[2] = rho[0]*(p[2] / p[0])**(1.0/gamma)
    rho[3] = rho[4] * (p[3] + Gamma*p[4]) / (p[4] + Gamma*p[3])
    
    n = [r/mass for r in rho]
    T = [p[i] / (n[i] * kb) for i in range(5)]
    
    # calc x
    x = np.array([0.0, args.L*0.5 - args.t*u[1], args.L*0.5, args.L*0.5 + args.t*u[3], args.L*0.5 + args.t*u[3] + args.t*u[4], args.L])
    
    if args.plt:
        plot_val(x, rho, "rho")
        plot_val(x, n, "n")
        plot_val(x, p, "p")
        plot_val(x, T, "T")
        
    if args.w:
        print("writing to file " + args.w + "_X.csv")
        write_values(x, rho, "rho", args.w)
        write_values(x, p, "p", args.w)
        write_values(x, n, "n", args.w)
        write_values(x, T, "T", args.w)
    
    print("done")

#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to perform test on the gradient function
    in numpy.

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep             = os.sep
    pwd             = os.getcwd()
    media_path      = pwd + "%cke_transport_terms%cmedia%c" %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    xf          = 2.0*pi
    xs          = 0.0
    k           = 1.0
    N           = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1048, 2096, 4192,\
                    8384, 16768]
    error2      = np.zeros(len(N))
    error1      = np.zeros(len(N))
    dx          = np.zeros(len(N))
    #---------------------------------------------------------------------#
    # Spatial step loop                                                   #
    #---------------------------------------------------------------------#
    for i, n in enumerate(N):
        print(n)
        dx[i]       = (xf-xs)/(n)
        x           = np.linspace(xs, xf, n+1)
        f           = np.cos(k*pi*x)
        exact       = -k*pi*np.sin(k*pi*x)
        approx2     = np.gradient(f, dx[i], edge_order=2)
        approx1     = np.gradient(f, dx[i], edge_order=1)
        error2[i]   = np.amax(abs(exact - approx2))
        error1[i]   = np.amax(abs(exact - approx1))
    #---------------------------------------------------------------------#
    # Plotting solution                                                   #
    #---------------------------------------------------------------------#
    plt.plot(x, exact, 'r', lw=1.5, label="Exact")
    plt.plot(x, approx2, 'bo', markevery=100, label="Approximate")
    plt.xlabel("$0 \leq x \leq 2.0\pi$")
    plt.ylabel("$\partial f / \partial x$")
    plt.grid(True)
    plt.legend(loc=0)
    plt.ylim([-4.5, 4.5])
    plt.savefig(media_path + "gradient-study-1.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting error                                                      #
    #---------------------------------------------------------------------#
    plt.loglog(dx, error2, 'r', lw=1.5, label="2nd order approximation")
    plt.loglog(dx, error1, 'b', lw=1.5, label="1st order approximation")
    plt.loglog(dx, 0.25*dx**2.0 , 'k', lw=1.5, label="$\sim c_{1} x^{2}$")
    plt.xlabel("Step size")
    plt.ylabel("Error")
    plt.grid(True)
    plt.legend(loc=0)
    plt.savefig(media_path + "gradient-study-2.pdf")
    plt.clf()

    print("**** Successful Run")
    sys.exit(0)

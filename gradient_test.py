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
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    N           = 1000
    xs          = 0.0
    xf          = 0.5*pi
    dx          = (xf-xs)/N
    x           = np.linspace(xs, xf, N+1)
    y           = np.linspace(xs, xf, N+1)
    [X1, X2]    = np.meshgrid(x, y)
    f           = np.zeros((N+1, N+1))
    dfE         = np.zeros((N+1, N+1))
    #---------------------------------------------------------------------#
    # Defining f and df/dx                                                #
    #---------------------------------------------------------------------#
    for j, yc in enumerate(y):
        for i, xc in enumerate(x):
            f[i,j]      = np.sin(pi*xc)*np.cos(pi*yc)
            dfE[i,j]    = -pi*np.sin(pi*xc)*np.sin(pi*yc)
    dfA = np.gradient(f, dx, edge_order=2)[1]
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, f, 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.show()
    #---------------------------------------------------------------------#
    # Plotting exact and approximate solutions                            #
    #---------------------------------------------------------------------#
    plt.subplot(2,1,1)
    cnt = plt.contourf(X1, X2, dfE, 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.title("Exact solution")
    plt.xlabel("$0 \leq x \leq \pi$") 
    plt.ylabel("$0 \leq y \leq \pi$") 
    plt.subplot(2,1,2)
    cnt = plt.contourf(X1, X2, dfA, 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.title("Approximate solution")
    plt.xlabel("$0 \leq x \leq \pi$") 
    plt.ylabel("$0 \leq y \leq \pi$") 
    plt.show()
    #---------------------------------------------------------------------#
    # Plotting error                                                      #
    #---------------------------------------------------------------------#
    err = abs(dfA - dfE)
    cnt = plt.contourf(X1, X2, err, 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.show()
    #---------------------------------------------------------------------#
    # Printing maximum error                                              #
    #---------------------------------------------------------------------#
    print("maximum error = %.8e"            %(np.amax(err)))

    print("**** Successful Run ****")
    sys.exit(0)

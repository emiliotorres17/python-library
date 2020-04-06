#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the strain rates of a
    fluid flied.

    **** Note:
                St[0]   --> S_{11}
                St[1]   --> S_{12}
                St[2]   --> S_{13}
                St[3]   --> S_{22}
                St[4]   --> S_{23}
                St[5]   --> S_{33}

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
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Strain rates                                                            #
#-------------------------------------------------------------------------#
def strain_rates(
        u1,                     # velocity-1 component
        u2,                     # velocity-2 component
        u3,                     # velocity-3 component
        h       = False,        # spatial step size
        dim     = False):       # number of spatial steps 

    """ Calculating the strain rates from the velocity fields """
    #---------------------------------------------------------------------#
    # Default settings                                                    #
    #---------------------------------------------------------------------#
    if h is False:
        pi  = np.pi
        N   = 64
        h   = (2.0*pi)/N
    if dim is False:
        dim = 64
    #---------------------------------------------------------------------#
    # Strain rates                                                        #
    #---------------------------------------------------------------------#
    St      = np.empty((6, dim, dim, dim))
    St[0]   =  np.gradient(u1, h, edge_order=2)[2]              # S11
    St[1]   = 0.5*(np.gradient(u1, h, edge_order=2)[1] +\
                   np.gradient(u2, h, edge_order=2)[2])         # S12
    St[2]   = 0.5*(np.gradient(u1, h, edge_order=2)[0] +\
                   np.gradient(u3, h, edge_order=2)[2])         # S13
    St[3]   =  np.gradient(u2, h, edge_order=2)[1]              # S22
    St[4]   = 0.5*(np.gradient(u2, h, edge_order=2)[0] +\
                   np.gradient(u3, h, edge_order=2)[1])         # S23
    St[5]   = (np.gradient(u3, h, edge_order=2)[0])             # S33

    return St
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(["clear"])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + "%cdata%c"          %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 256
    x0          = 0.0
    xf          = 1.0
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    dx          = (xf-x0)/N
    [X1, X2]    = np.meshgrid(x,y)
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating velocities                                            #
    #---------------------------------------------------------------------#
    ux          = np.zeros((N+1, N+1, N+1))
    uy          = np.zeros((N+1, N+1, N+1))
    uz          = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the velocities                                          #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                ux[i,j,k]   = x[i]*y[j]*z[k]
                uy[i,j,k]   = x[i] + 3.0*y[j] + z[k]
                uz[i,j,k]   = x[i]*y[j]*z[k]**2.0
        #-----------------------------------------------------------------#
        # Print count                                                     #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate solution                                #
    #---------------------------------------------------------------------#
    approx  = strain_rates(ux, uy, uz, dx, N+1)[0]
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del ux
    del uy
    del uz
    #=====================================================================#
    # Exact solution                                                      #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating strain rates                                          #
    #---------------------------------------------------------------------#
    s11     = np.zeros((N+1, N+1, N+1))
    s12     = np.zeros((N+1, N+1, N+1))
    s13     = np.zeros((N+1, N+1, N+1))
    s22     = np.zeros((N+1, N+1, N+1))
    s23     = np.zeros((N+1, N+1, N+1))
    s33     = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the exact strain rates                                  #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                s11[i,j,k]  = y[j]*z[k]
                s12[i,j,k]  = 0.5*(x[i]*z[k] + 1.0)
                s13[i,j,k]  = 0.5*(x[i]*y[j] + y[j]*z[k]**2.0)
                s22[i,j,k]  = 3.0
                s23[i,j,k]  = 0.5*(x[i]*z[k]**2.0 + 1.0)
                s33[i,j,k]  = 2.0*x[i]*y[j]*z[k]
        #-----------------------------------------------------------------#
        # Print count                                                     #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #=====================================================================#
    # Error                                                               #
    #=====================================================================#
    error   = abs(approx - s11)
    #=====================================================================#
    # Plotting solutions                                                  #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, s11[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"    %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"    %(x0, xf))
    plt.colorbar()
    plt.show()
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, approx[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"    %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"    %(x0, xf))
    plt.colorbar()
    plt.show()
    #---------------------------------------------------------------------#
    # Error                                                               #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"    %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"    %(x0, xf))
    plt.colorbar()
    plt.show()

    print("**** Successful Run ****")
    sys.exit(0)

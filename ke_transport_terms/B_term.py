#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate the B-term in the kinetic
    energy transport equation

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
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from strain_rates   import strain_rates
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# B term                                                                  #
#-------------------------------------------------------------------------#
def b_term(
        U1,                 # velocity-1 component
        U2,                 # velocity-2 component
        U3,                 # velocity-3 component
        h,                  # spatial step size
        NU):                # viscosity

    """ Calculating the B term in the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33)  = strain_rates(U1, U2, U3, h)
    #---------------------------------------------------------------------#
    # Calculating the 9 gradient terms                                    #
    #---------------------------------------------------------------------#
    Term1   = np.gradient(U1*S11, h, edge_order=2)[0]
    Term2   = np.gradient(U1*S12, h, edge_order=2)[1]
    Term3   = np.gradient(U1*S13, h, edge_order=2)[2]
    Term4   = np.gradient(U2*S12, h, edge_order=2)[0]
    Term5   = np.gradient(U2*S22, h, edge_order=2)[1]
    Term6   = np.gradient(U2*S23, h, edge_order=2)[2]
    Term7   = np.gradient(U3*S13, h, edge_order=2)[0]
    Term8   = np.gradient(U3*S23, h, edge_order=2)[1]
    Term9   = np.gradient(U3*S33, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating B term                                                  #
    #---------------------------------------------------------------------#
    B       = 2.0*NU*(Term1 + Term2 + Term3 + Term4 + Term5 + Term6 +\
                        Term7 + Term8 + Term9)
    #---------------------------------------------------------------------#
    # Deleting term variables                                             #
    #---------------------------------------------------------------------#
    del Term1
    del Term2
    del Term3
    del Term4
    del Term5
    del Term6
    del Term7
    del Term8
    del Term9

    return B
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep         = os.sep
    pwd         = os.getcwd()
    media_path  = pwd + "%cmedia%c"         %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 400
    pi          = np.pi
    x0          = 0.0
    xf          = 0.5
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    dx          = (xf - x0)/N
    nu          = 1.0
    #=====================================================================#
    # Preallocating space                                                 #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Velocities                                                          #
    #---------------------------------------------------------------------#
    ux          = np.zeros((N+1, N+1, N+1))
    uy          = np.zeros((N+1, N+1, N+1))
    uz          = np.zeros((N+1, N+1, N+1))
    #=====================================================================#
    # Defining velocities and strain rates                                #
    #=====================================================================#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                #---------------------------------------------------------#
                # Velocities                                              #
                #---------------------------------------------------------#
                ux[i,j,k]   = np.cos(pi*x[i])
                uy[i,j,k]   = np.sin(2.0*pi*y[j])
                uz[i,j,k]   = np.sin(pi*z[k])
        #-----------------------------------------------------------------#
        # print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #=====================================================================#
    # Calculating the solution                                            #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Approximation solution                                              #
    #---------------------------------------------------------------------#
    B_approx    = b_term(ux, uy, uz, dx, nu)
    del ux
    del uy
    del uz
    del dx
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    term1       = np.zeros((N+1, N+1, N+1))
    term5       = np.zeros((N+1, N+1, N+1))
    term9       = np.zeros((N+1, N+1, N+1))
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                #---------------------------------------------------------#
                # Exact solution terms                                    #
                #---------------------------------------------------------#
                term1[i,j,k]    = pi**2.0 - 2.0*pi**2.0*(np.cos(pi*x[i]))**2.0
                term5[i,j,k]    = 8.0*pi**2.0*(np.cos(2.0*pi*y[j]))**2.0 - 4.0*pi**2.0
                term9[i,j,k]    = 2.0*pi**2.0*(np.cos(pi*z[k]))**2.0 - pi**2.0

        #-----------------------------------------------------------------#
        # print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    B_exact     = 2.0*nu*(term1 + term5 + term9)
    #---------------------------------------------------------------------#
    # Calculating the error                                               #
    #---------------------------------------------------------------------#
    error   = abs(B_exact - B_approx)
    #=====================================================================#
    # Plotting the results                                                #
    #=====================================================================#
    [X1, X2]    = np.meshgrid(x, y)
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, B_exact[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$-0.5 \leq x_{1} \leq 0.5$")
    plt.ylabel("$-0.5 \leq x_{2} \leq 0.5$")
    plt.colorbar()
    plt.savefig(media_path + "B-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximation solution                                              #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, B_approx[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$-0.5 \leq x_{1} \leq 0.5$")
    plt.ylabel("$-0.5 \leq x_{2} \leq 0.5$")
    plt.colorbar()
    plt.savefig(media_path + "B-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$-0.5 \leq x_{1} \leq 0.5$")
    plt.ylabel("$-0.5 \leq x_{2} \leq 0.5$")
    plt.colorbar()
    plt.savefig(media_path + "B-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

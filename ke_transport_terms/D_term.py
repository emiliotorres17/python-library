#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the D term in the
    kinetic energy transport equation.

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
# D term                                                                  #
#-------------------------------------------------------------------------#
def d_term(
        U1,                 # velocity-1 component
        U2,                 # velocity-2 component
        U3,                 # velocity-3 component
        NU,                 # viscosity
        h):                 # spatial step size

    """ Calculating the D term in the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33)  = strain_rates(U1, U2, U3, h)
    #---------------------------------------------------------------------#
    # Calculating the 6 terms                                             #
    #---------------------------------------------------------------------#
    Term1       = np.multiply(S11, S11)
    Term2       = 2.0*np.multiply(S12, S12)
    Term3       = 2.0*np.multiply(S13, S13)
    Term4       = np.multiply(S22, S22)
    Term5       = 2.0*np.multiply(S23, S23)
    Term6       = np.multiply(S33, S33)
    #---------------------------------------------------------------------#
    # Calculating the D term                                              #
    #---------------------------------------------------------------------#
    D           = -2.0*NU*(Term1 + Term2 + Term3 + Term4 + Term5 + Term6)

    return D
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
    media_path  = pwd + "%cmedia%c"         %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 128
    x0          = 0.0
    xf          = 1.0
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x, y)
    dx          = (xf - x0)/N
    nu          = 1.0
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating velocities                                            #
    #---------------------------------------------------------------------#
    ux  = np.zeros((N+1, N+1, N+1))
    uy  = np.zeros((N+1, N+1, N+1))
    uz  = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the velocities                                          #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                ux[i,j,k]   = x[i]**2.0*y[j]*z[k]
                uy[i,j,k]   = x[i]*y[j]**2.0*z[k]
                uz[i,j,k]   = x[i]*y[j]*z[k]**2.0
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate solution                                #
    #---------------------------------------------------------------------#
    D_approx    = d_term(ux, uy, uz, nu, dx)
    print("***** Minimum value")
    print(np.amax(ux))
    print(np.amax(uy))
    print(np.amax(uz))
    #---------------------------------------------------------------------#
    # Clearing variables
    #---------------------------------------------------------------------#
    del ux
    del uy
    del uz
    #=====================================================================#
    # Exact solution                                                      #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating terms                                                 #
    #---------------------------------------------------------------------#
    term1       = np.zeros((N+1, N+1, N+1))
    term2       = np.zeros((N+1, N+1, N+1))
    term3       = np.zeros((N+1, N+1, N+1))
    term4       = np.zeros((N+1, N+1, N+1))
    term5       = np.zeros((N+1, N+1, N+1))
    term6       = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the terms in the D term                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                term1[i,j,k]    = 4.0*x[i]**2.0*y[j]**2.0*z[k]**2.0
                term2[i,j,k]    = (0.5*(x[i]**2.0*z[k]+y[j]**2.0*z[k]))**2.0
                term3[i,j,k]    = (0.5*(x[i]**2.0*y[j]+y[j]*z[k]**2.0))**2.0
                term4[i,j,k]    = 4.0*x[i]**2.0*y[j]**2.0*z[k]**2.0
                term5[i,j,k]    = (0.5*(x[i]*y[j]**2.0+x[i]*z[k]**2.0))**2.0
                term6[i,j,k]    = 4.0*x[i]**2.0*y[j]**2.0*z[k]**2.0
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the D term                                              #
    #---------------------------------------------------------------------#
    D_exact     = -2.0*nu*(term1 + 2.0*term2 + 2.0*term3 + term4 +\
                            2.0*term5 + term6)
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del term1
    del term2
    del term3
    del term4
    del term5
    del term6
    #=====================================================================#
    # Error                                                               #
    #=====================================================================#
    error   = abs(D_exact - D_approx)
    #=====================================================================#
    # Plotting the solutions                                              #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, D_exact[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "D-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, D_approx[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "D-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error                                                               #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "D-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

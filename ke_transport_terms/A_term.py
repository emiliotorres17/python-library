#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the A term in the
    kinetic energy transport term.

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
# Terms in the A term                                                     #
#-------------------------------------------------------------------------#
def a_terms(
        U1,                 # velocity 1 field
        U2,                 # velocity 2 field
        U3,                 # velocity 3 field
        P,                  # pressure field
        h,                  # spatial step size
        rho):               # density

    """ Calculating the terms in the A term of the kinetic energy transport
    equation """
    #---------------------------------------------------------------------#
    # Calculating the 3 different terms                                   #
    #---------------------------------------------------------------------#
    term1       = (1.0/rho)*np.gradient(P*U1, h, edge_order=2)[0]
    term2       = (1.0/rho)*np.gradient(P*U2, h, edge_order=2)[1]
    term3       = (1.0/rho)*np.gradient(P*U3, h, edge_order=2)[2]

    return term1, term2, term3
#-------------------------------------------------------------------------#
# A term                                                                  #
#-------------------------------------------------------------------------#
def a_term(
        U1,                 # velocity-1 field
        U2,                 # velocity-2 field
        U3,                 # velocity-3 field
        P,                  # pressure field
        den,                # density field
        h):                 # step size

    """ Calculating the A term in the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Calculating the 3 different terms                                   #
    #---------------------------------------------------------------------#
    term1       = np.gradient(P*U1, h, edge_order=2)[0]
    term2       = np.gradient(P*U2, h, edge_order=2)[1]
    term3       = np.gradient(P*U3, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating the A term                                              #
    #---------------------------------------------------------------------#
    a   = -(1.0/den)*(term1 + term2 + term3)

    return a
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Preamble                                                            #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep         = os.sep
    pi          = np.pi
    #---------------------------------------------------------------------#
    # Defining paths                                                      #
    #---------------------------------------------------------------------#
    media_path      = "media%c"                     %(sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 400
    x0          = -0.5
    xf          = 0.5
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x, y)
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    p           = np.zeros((N+1, N+1, N+1))
    dpdx        = np.zeros((N+1, N+1, N+1))
    dpdy        = np.zeros((N+1, N+1, N+1))
    dpdz        = np.zeros((N+1, N+1, N+1))
    ux          = np.zeros((N+1, N+1, N+1))
    duxdx       = np.zeros((N+1, N+1, N+1))
    uy          = np.zeros((N+1, N+1, N+1))
    duydy       = np.zeros((N+1, N+1, N+1))
    uz          = np.zeros((N+1, N+1, N+1))
    duzdz       = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                p[i,j,k]        = np.sin(2.0*pi*x[i])*np.cos(3.0*pi*y[j])*(np.sin(2.0*pi*z[k]))**2.0
                ux[i,j,k]       = np.cos(pi*x[i])*np.sin(3.0*pi*z[k])
                uy[i,j,k]       = np.sin(pi*y[j])
                uz[i,j,k]       = np.cos(3.0*pi*y[j])*np.cos(pi*z[k])
                dpdx[i,j,k]     = 2.0*pi*np.cos(2.0*pi*x[i])*np.cos(3.0*pi*y[j])*(np.sin(2.0*pi*z[k]))**2.0
                dpdy[i,j,k]     = -3.0*pi*np.sin(2.0*pi*x[i])*np.sin(3.0*pi*y[j])*(np.sin(2.0*pi*z[k]))**2.0
                dpdz[i,j,k]     = 4.0*pi*np.sin(2.0*pi*x[i])*np.cos(3.0*pi*y[j])*np.sin(2.0*pi*z[k])*np.cos(2.0*pi*z[k])
                duxdx[i,j,k]    = -pi*np.sin(pi*x[i])*np.sin(3.0*pi*z[k])
                duydy[i,j,k]    = pi*np.cos(pi*y[j])
                duzdz[i,j,k]    = -pi*np.cos(3.0*pi*y[j])*np.sin(pi*z[k])
        if print_count == 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the exact A term                                        #
    #---------------------------------------------------------------------#
    A_exact     = -1.0*(p*duxdx + ux*dpdx + p*duydy + uy*dpdy + p*duzdz + uz*dpdz)
    #---------------------------------------------------------------------#
    # Calculating the approximate A term                                  #
    #---------------------------------------------------------------------#
    A_approx    = a_term(ux, uy, uz, p, 1.0, dx)
    #---------------------------------------------------------------------#
    # Calculating the error                                               #
    #---------------------------------------------------------------------#
    error       = abs(A_approx - A_exact)
    #---------------------------------------------------------------------#
    # Plotting A exact solution                                           #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, A_exact[:,:,31], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting A approximate solution                                     #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, A_approx[:,:,31], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting the error                                                  #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,31], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

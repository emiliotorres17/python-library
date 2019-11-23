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
def a_term(
        U1,                 # velocity-1 field
        U2,                 # velocity-2 field
        U3,                 # velocity-3 field
        p,                  # pressure field
        rho,                # density field
        h):                 # step size

    """ Calculating the A term in the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Calculating the 3 different terms                                   #
    #---------------------------------------------------------------------#
    term1       = np.gradient(p*U1, h, edge_order=2)[0]
    term2       = np.gradient(p*U2, h, edge_order=2)[1]
    term3       = np.gradient(p*U3, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating the A term                                              #
    #---------------------------------------------------------------------#
    a   = -(1.0/rho)*(term1 + term2 + term3)

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
    data_path       = "data%c"                      %(sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    dx          = (2.0*pi)/(64.0)
    x1          = np.linspace(0.0, 2.0*pi, 64)
    x2          = np.linspace(0.0, 2.0*pi, 64)
    [X1, X2]    = np.meshgrid(x1, x2)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    u1          = np.load(data_path + "velocity1.npy")
    u2          = np.load(data_path + "velocity2.npy")
    u3          = np.load(data_path + "velocity3.npy")
    pressure    = np.load(data_path + "pressure.npy")
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    u1          = u1[:,:,:, 5]
    u2          = u2[:,:,:, 5]
    u3          = u3[:,:,:, 5]
    pressure    = pressure[:,:,:, 5]
    #---------------------------------------------------------------------#
    # Calculating the A-term                                              #
    #---------------------------------------------------------------------#
    A   = a_term(u1, u2, u3, pressure, 1.0, dx)
    #---------------------------------------------------------------------#
    # Generating the contour plot                                         #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, A[:,:,32], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel(r"$0 \leq x_{1} \leq 2\pi$")
    plt.ylabel(r"$0 \leq x_{2} \leq 2\pi$")
    plt.colorbar()
    plt.show()

    print("**** Successful Run ****")
    sys.exit(0)

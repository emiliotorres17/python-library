#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the strain rates of a
    fluid flied.

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
import scipy.io as sci
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
        h):                     # spatial step size

    """ Calculating the strain rates from the velocity fields """
    #---------------------------------------------------------------------#
    # Strain rates                                                        #
    #---------------------------------------------------------------------#
    s11     =  np.gradient(u1, h, edge_order=2)[0]              # s11
    s12     = 0.5*(np.gradient(u1, h, edge_order=2)[1] +\
                    np.gradient(u2, h, edge_order=2)[0])        # s12
    s13     = 0.5*(np.gradient(u1, h, edge_order=2)[2] +\
                    np.gradient(u3, h, edge_order=2)[0])        # s13
    s22     =  np.gradient(u2, h, edge_order=2)[1]              # s22
    s23     = 0.5*(np.gradient(u2, h, edge_order=2)[2] +\
                    np.gradient(u3, h, edge_order=2)[1])        # s23
    s33     = (np.gradient(u3, h, edge_order=2)[2])             # s33

    return s11, s12, s13, s22, s23, s33
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
    data_path   = pwd + "%cdata%c"          %(sep, sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    vel         = sci.loadmat(data_path + 'velocities-20.mat')
    vel1        = vel['u1']
    vel2        = vel['u2']
    vel3        = vel['u3']
    print("**** Loaded data ****")
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    pi          = np.pi
    dx          = (2.0*pi)/64.0
    x1          = np.linspace(0, 2.0*pi, 64)
    x2          = np.linspace(0, 2.0*pi, 64)
    [X1, X2]    = np.meshgrid(x1, x2)
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33) = strain_rates(vel1, vel2, vel3, dx)
    #---------------------------------------------------------------------#
    # Making contours                                                     #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, S23[:,:,44], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel(r"$0 \leq x_{1} \leq 2\pi$")
    plt.ylabel(r"$0 \leq x_{2} \leq 2\pi$")
    plt.colorbar()
    plt.show()
    for i in range(0, 5):
        print(S23[i,12,45])

    print("**** Successful Run ****")
    sys.exit(0)

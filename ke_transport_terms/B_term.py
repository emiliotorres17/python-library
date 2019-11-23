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
        u1,                 # velocity-1 component
        u2,                 # velocity-2 component
        u3,                 # velocity-3 component
        s11,                # strain rate-11 component
        s12,                # strain rate-12 component
        s13,                # strain rate-13 component
        s22,                # strain rate-22 component
        s23,                # strain rate-23 component
        s33,                # strain rate-33 component 
        h,                  # spatial step size
        nu):                # viscosity

    """ Calculating the B term in the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Finding the 9 gradient terms                                        #
    #---------------------------------------------------------------------#
    term1   = np.gradient(u1*s11, dx, edge_order=2)[0]
    term2   = np.gradient(u1*s12, dx, edge_order=2)[1]
    term3   = np.gradient(u1*s13, dx, edge_order=2)[2]
    term4   = np.gradient(u2*s12, dx, edge_order=2)[0]
    term5   = np.gradient(u2*s22, dx, edge_order=2)[1]
    term6   = np.gradient(u2*s23, dx, edge_order=2)[2]
    term7   = np.gradient(u3*s13, dx, edge_order=2)[0]
    term8   = np.gradient(u3*s23, dx, edge_order=2)[1]
    term9   = np.gradient(u3*s33, dx, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating B term                                                  #
    #---------------------------------------------------------------------#
    b       = 2.0*nu*(term1 + term2 + term3 + term4 + term5 + term6 +\
                        term7 + term8 + term9)

    return b
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
    data_path   = pwd + "%cdata%c"                  %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    pi          = np.pi
    x1          = np.linspace(0, 2.0*pi, 64)
    x2          = np.linspace(0, 2.0*pi, 64)
    [X1, X2]    = np.meshgrid(x1, x2)
    dx          = (2.0*pi)/64.0
    nu          = 0.000185
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    u1          = np.load(data_path + 'velocity1.npy')
    u2          = np.load(data_path + 'velocity2.npy')
    u3          = np.load(data_path + 'velocity3.npy')
    #---------------------------------------------------------------------#
    # Defining velocities and strain rates                                #
    #---------------------------------------------------------------------#
    u1          = u1[:,:,:,50]
    u2          = u2[:,:,:,50]
    u3          = u3[:,:,:,50]
    (s11, s12, s13, s22, s23, s33)  = strain_rates(u1, u2, u3, dx) 
    #---------------------------------------------------------------------#
    # Calculating the B term                                              #
    #---------------------------------------------------------------------#
    b           = b_term(u1, u2, u3, s11, s12, s13, s22, s23, s33, dx, nu)
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    cnt         = plt.contourf(X1, X2, b[:,:,32], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.show()

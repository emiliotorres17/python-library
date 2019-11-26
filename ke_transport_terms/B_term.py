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

#
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
    term1   = np.gradient(U1*S11, h, edge_order=2)[0]
    term2   = np.gradient(U1*S12, h, edge_order=2)[1]
    term3   = np.gradient(U1*S13, h, edge_order=2)[2]
    term4   = np.gradient(U2*S12, h, edge_order=2)[0]
    term5   = np.gradient(U2*S22, h, edge_order=2)[1]
    term6   = np.gradient(U2*S23, h, edge_order=2)[2]
    term7   = np.gradient(U3*S13, h, edge_order=2)[0]
    term8   = np.gradient(U3*S23, h, edge_order=2)[1]
    term9   = np.gradient(U3*S33, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating B term                                                  #
    #---------------------------------------------------------------------#
    B       = 2.0*NU*(term1 + term2 + term3 + term4 + term5 + term6 +\
                        term7 + term8 + term9)

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
    u1          = u1[:,:,:,100]
    u2          = u2[:,:,:,100]
    u3          = u3[:,:,:,100]
    (s11, s12, s13, s22, s23, s33)  = strain_rates(u1, u2, u3, dx)
    #---------------------------------------------------------------------#
    # Calculating the B term                                              #
    #---------------------------------------------------------------------#
    b           = b_term(u1, u2, u3,  dx, nu)
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    cnt         = plt.contourf(X1, X2, b[:,:,32], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.show()

    print("**** Successful Run ****")
    sys.exit(0)

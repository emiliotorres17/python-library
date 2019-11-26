#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the P term in the
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
# P term                                                                  #
#-------------------------------------------------------------------------#
def p_term(
        U1,                     # velocity-1 component
        U2,                     # velocity-2 component
        U3,                     # velocity-3 component
        Tau11,                  # tau-11 component
        Tau12,                  # tau-12 component
        Tau13,                  # tau-13 component
        Tau22,                  # tau-22 component
        Tau23,                  # tau-23 component
        Tau33,                  # tau-33 component
        h):                     # spatial step size

    """ Calculating the P term """
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33)  = strain_rates(U1, U2, U3, h)
    #---------------------------------------------------------------------#
    # Calculating the terms in the P term                                 #
    #---------------------------------------------------------------------#
    term1   = np.multiply(Tau11, S11)
    term2   = 2.0*np.multiply(Tau12, S12)
    term3   = 2.0*np.multiply(Tau13, S13)
    term4   = np.multiply(Tau22, S22)
    term5   = 2.0*np.multiply(Tau23, S23)
    term6   = np.multiply(Tau33, S33)
    #---------------------------------------------------------------------#
    # Calculating the P term                                              #
    #---------------------------------------------------------------------#
    P   = term1 + term2 + term3 + term4 + term5 + term6

    return P
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
    data_path   = pwd + "%cdata%c"              %(sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    x1          = np.linspace(0.0, 2.0*pi, 64)
    x2          = np.linspace(0.0, 2.0*pi, 64)
    [X1, X2]    = np.meshgrid(x1, x2)
    dx          = (2.0*pi)/64.0
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    u1      = np.load(data_path + "velocity1.npy")
    u2      = np.load(data_path + "velocity2.npy")
    u3      = np.load(data_path + "velocity3.npy")
    tau11   = np.load(data_path + "tau11.npy")
    tau12   = np.load(data_path + "tau12.npy")
    tau13   = np.load(data_path + "tau13.npy")
    tau22   = np.load(data_path + "tau22.npy")
    tau23   = np.load(data_path + "tau23.npy")
    tau33   = np.load(data_path + "tau33.npy")
    #---------------------------------------------------------------------#
    # Defining specific time step data                                    #
    #---------------------------------------------------------------------#
    tslice  = 25
    u1      = u1[:,:,:,tslice]
    u2      = u2[:,:,:,tslice]
    u3      = u3[:,:,:,tslice]
    tau11   = tau11[:,:,:,tslice]
    tau12   = tau12[:,:,:,tslice]
    tau13   = tau13[:,:,:,tslice]
    tau22   = tau22[:,:,:,tslice]
    tau23   = tau23[:,:,:,tslice]
    tau33   = tau33[:,:,:,tslice]
    #---------------------------------------------------------------------#
    # Calculating the P term                                              #
    #---------------------------------------------------------------------#
    p   = p_term(u1, u2, u3, tau11, tau12, tau13, tau22, tau23, tau33, dx)
    #---------------------------------------------------------------------#
    # Plotting P results                                                  #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, p[:,:,32], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("0 \leq x_{1} \leq 2 \pi")
    plt.ylabel("0 \leq x_{2} \leq 2 \pi")
    plt.colorbar()
    plt.show()

    print("**** Successful Run ****")
    sys.exit(0)

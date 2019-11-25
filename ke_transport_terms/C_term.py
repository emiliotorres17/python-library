#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the C term in the
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
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# C term                                                                  #
#-------------------------------------------------------------------------#
def c_term(
        U1,                 # velocity-1 component
        U2,                 # velocity-2 component
        U3,                 # velocity-3 component
        Tau11,              # tau-11 component
        Tau12,              # tau-12 component
        Tau13,              # tau-13 component
        Tau22,              # tau-22 component
        Tau23,              # tau-23 component
        Tau33,              # tau-33 component
        h):                 # spatial step size

    """ Calculating the C term """
    #---------------------------------------------------------------------#
    # Calculating the terms in the C term                                 #
    #---------------------------------------------------------------------#
    term1   = np.gradient(Tau11*U1, h, edge_order=2)[0]
    term2   = np.gradient(Tau12*U1, h, edge_order=2)[1]
    term3   = np.gradient(Tau13*U1, h, edge_order=2)[2]
    term4   = np.gradient(Tau12*U2, h, edge_order=2)[0]
    term5   = np.gradient(Tau22*U2, h, edge_order=2)[1]
    term6   = np.gradient(Tau23*U2, h, edge_order=2)[2]
    term7   = np.gradient(Tau13*U3, h, edge_order=2)[0]
    term8   = np.gradient(Tau23*U3, h, edge_order=2)[1]
    term9   = np.gradient(Tau33*U3, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating the C term                                              #
    #---------------------------------------------------------------------#
    C   = -(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8\
                + term9)

    return C
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
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    pi          = np.pi
    x1          = np.linspace(0, 2.0*pi, 64)
    x2          = np.linspace(0, 2.0*pi, 64)
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
    # Defining plotting variables                                         #
    #---------------------------------------------------------------------#
    tslice  = 20
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
    # Calculating the c term                                              #
    #---------------------------------------------------------------------#
    c       = c_term(u1, u2, u3, tau11, tau12, tau13, tau22, tau23, tau33,\
                        dx)
    #---------------------------------------------------------------------#
    # Plotting the results                                                #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, c[:,:,32], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 2\pi$")
    plt.ylabel("$0 \leq x_{1} \leq 2\pi$")
    plt.colorbar()
    plt.show()

    print("**** Successful Run ****")
    sys.exit(0)

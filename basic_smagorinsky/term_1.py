#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate term 1 in the enstrophy
    transport study.

    (nu_{sgs}/nu)(B_{Omega} + D_{Omega})

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import sys
from subprocess import call
import numpy as np
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Term 1                                                                  #
#-------------------------------------------------------------------------#
def term_1(
        omega1,                 # vorticity-1
        omega2,                 # vorticity-2
        omega3,                 # vorticity-3
        enst,                   # enstrophy
        nu_sgs,                 # turbulent viscosity
        h       = True):        # spatial step size

    """ Calculating term 1 in the SGS enstrophy study """
    #---------------------------------------------------------------------#
    # Setting default values                                              #
    #---------------------------------------------------------------------#
    if h is True:
        h = 2.0*np.pi/64.0
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    term    = np.zeros((64,64,64))
    #---------------------------------------------------------------------#
    # Enstrophy term                                                      #
    #---------------------------------------------------------------------#
    term    += np.gradient(\
                np.gradient(enst, h, edge_order=2)[2], h, edge_order=2)[2]
    term    += np.gradient(\
                np.gradient(enst, h, edge_order=2)[1], h, edge_order=2)[1]
    term    += np.gradient(\
                np.gradient(enst, h, edge_order=2)[0], h, edge_order=2)[0]
    #---------------------------------------------------------------------#
    # Dissipation                                                         #
    #---------------------------------------------------------------------#
    omega1_grad = np.gradient(omega1, h, edge_order=2)
    omega2_grad = np.gradient(omega2, h, edge_order=2)
    omega3_grad = np.gradient(omega3, h, edge_order=2)
    term        -= np.square(omega1_grad[2])
    term        -= np.square(omega1_grad[1])
    term        -= np.square(omega1_grad[0])
    term        -= np.square(omega2_grad[2])
    term        -= np.square(omega2_grad[1])
    term        -= np.square(omega2_grad[0])
    term        -= np.square(omega3_grad[2])
    term        -= np.square(omega3_grad[1])
    term        -= np.square(omega3_grad[0])
    #---------------------------------------------------------------------#
    # Applying the subgrid stress                                         #
    #---------------------------------------------------------------------#
    term        *= nu_sgs

    return term
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])

    print('**** This has not been unit test *****')

    sys.exit(0)

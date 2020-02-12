#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the SGS transport in the
    enstrophy transport equation.

    **** Note: The spectral differentiation has not been added yet.

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
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from enstrophy_transport_term.psi_enstrophy     import psi
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
def SGS_transport(
        w1,                 # vorticity-1 component
        w2,                 # vorticity-2 component
        w3,                 # vorticity-3 component
        Tau,                # SGS; (6,64,64,64)
        h       = False,    # spatial step size
        flag    = True):    # spectral flag; default is gradient tool

    """ Calculating the SGS transport term """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if h is False:
        Pi  = np.pi
        num = 64
        h   = (2.0*Pi)/num
    #---------------------------------------------------------------------#
    # Preallocating  variables                                            #
    #---------------------------------------------------------------------#
    dim         = np.shape(w1)[0]
    SGS_trans   = np.zeros(dim, dim, dim)
    #---------------------------------------------------------------------#
    # Calculating Psi                                                     #
    #---------------------------------------------------------------------#
    Psi = psi(Tau, h, True)
    #---------------------------------------------------------------------#
    # Calculating SGS transport                                           #
    #---------------------------------------------------------------------#
    if flag is False:       # spectral flag
        print("**** Error:Spectral differentiation has not been setup ****")
        sys.exit(1)
    else:                   # gradient tool flag
        #-----------------------------------------------------------------#
        # Terms 1-3 (i = 1)                                               #
        #-----------------------------------------------------------------#
        SGS_trans   += np.gradient(w1*Psi[0], h, edge_orer=2)[0]
        SGS_trans   += np.gradient(w1*Psi[1], h, edge_orer=2)[1]
        SGS_trans   += np.gradient(w1*Psi[2], h, edge_orer=2)[2]
        #-----------------------------------------------------------------#
        # Terms 4-6 (i = 2)                                               #
        #-----------------------------------------------------------------#
        SGS_trans   += np.gradient(w2*Psi[3], h, edge_orer=2)[0]
        SGS_trans   += np.gradient(w2*Psi[4], h, edge_orer=2)[1]
        SGS_trans   += np.gradient(w2*Psi[5], h, edge_orer=2)[2]
        #-----------------------------------------------------------------#
        # Terms 7-9 (i = 3)                                               #
        #-----------------------------------------------------------------#
        SGS_trans   += np.gradient(w3*Psi[6], h, edge_orer=2)[0]
        SGS_trans   += np.gradient(w3*Psi[7], h, edge_orer=2)[1]
        SGS_trans   += np.gradient(w3*Psi[8], h, edge_orer=2)[2]

    SGS_trans *= -1

    return SGS_trans
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    print("**** Has not been unit tested ****")

    sys.exit(0)

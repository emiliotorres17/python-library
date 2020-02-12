#!/usr/bin/env python
"""========================================================================
Purpose:
    The purpose of this script to calculate the C term, laplacian(Omega),
    in the enstrophy transport equation.

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
import matplotlob.pyplot as plt
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# C term                                                                  #
#-------------------------------------------------------------------------#
def c_term(
        enst,                           # enstrophy field
        dx      = False,                # spatial step size
        Nu      = False,                # viscosity
        flag    = True):                # spectral flag; default is gradient

    """ Calculating the diffusion term in enstrophy transport equation """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if dx is False:
        num = 64
        Pi  = np.pi
        dx  = (2.0*pi)/num
    if Nu is False:
        Nu  = 0.000185
    #---------------------------------------------------------------------#
    # Calculating the Laplacian                                           #
    #---------------------------------------------------------------------#
    diim    = enst.shape[0]
    if flag is True:
        term1   = np.gradient(np.gradient(enst, dx, edge_order=2)[0],\
                                dx, edge_order=2)[0]
        term2   = np.gradient(np.gradient(enst, dx, edge_order=2)[1],\
                                dx, edge_order=2)[1]
        term3   = np.gradient(np.gradient(enst, dx, edge_order=2)[2],\
                                dx, edge_order=2)[2]
    else:
        term1   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term1   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term1)).real
        term2   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term2   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term2)).real
        term3   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term3   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term3)).real
    #---------------------------------------------------------------------#
    # Calculating dissipation term                                        #
    #---------------------------------------------------------------------#
    c   = np.zeros((dim,dim,dim))
    c   += Nu*(term1 + term2 + term3)

    return c
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

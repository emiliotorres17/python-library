#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the Psi operator in the
    enstrophy transport equation.

    **** Note: This subroutine has not applied spectral differentiation

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
# Psi operator                                                            #
#-------------------------------------------------------------------------#
def psi_enstrophy(
        Tau,                    # SGS; (6,64,64,64)
        h       = False,        # spatial step size
        flag    = True):        # spectral flag; default is gradient tool

    """ Calculating the psi operator for the transport and production of
    the enstrophy """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if h is False:
        Pi  = np.pi
        N   = 64
        h   = (2.0*Pi)/N
    #---------------------------------------------------------------------#
    # Preallocation variables                                             #
    #---------------------------------------------------------------------#
    dim = np.shape(Tau)[1]
    Psi = np.zeros((9, dim, dim, dim))
    #---------------------------------------------------------------------#
    # Calculating psi using spectral methods                              #
    #---------------------------------------------------------------------#
    if flag is False:
        kspec   = np.fft.fftfreq(dim) * dim
        Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
        #-----------------------------------------------------------------#
        # Psi_{11}                                                        #
        #-----------------------------------------------------------------#
        Psi[0]  = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[2])).real -\
                    np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[1])).real
        #-----------------------------------------------------------------#
        # Psi_{12}                                                        #
        #-----------------------------------------------------------------#
        Psi[1]  = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[4])).real -\
                    np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[3])).real
        #-----------------------------------------------------------------#
        # Psi_{13}                                                        #
        #-----------------------------------------------------------------#
        Psi[2]  = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[5])).real -\
                    np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[4])).real
        #-----------------------------------------------------------------#
        # Psi_{21}                                                        #
        #-----------------------------------------------------------------#
        Psi[3]  = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[0])).real -\
                    np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[2])).real
        #-----------------------------------------------------------------#
        # Psi_{22}                                                        #
        #-----------------------------------------------------------------#
        Psi[4]  = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[1])).real -\
                    np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[4])).real
        #-----------------------------------------------------------------#
        # Psi_{23}                                                        #
        #-----------------------------------------------------------------#
        Psi[5]  = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(Tau[2])).real -\
                    np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[5])).real
        #-----------------------------------------------------------------#
        # Psi_{31}                                                        #
        #-----------------------------------------------------------------#
        Psi[6]  = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[1])).real -\
                    np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[0])).real
        #-----------------------------------------------------------------#
        # Psi_{32}                                                        #
        #-----------------------------------------------------------------#
        Psi[7]  = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[3])).real -\
                    np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[1])).real
        #-----------------------------------------------------------------#
        # Psi_{33}                                                        #
        #-----------------------------------------------------------------#
        Psi[8]  = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(Tau[4])).real -\
                    np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(Tau[2])).real
    #---------------------------------------------------------------------#
    # Calculating psi using gradient tool                                 #
    #---------------------------------------------------------------------#
    else:
        #-----------------------------------------------------------------#
        # Psi_{11}                                                        #
        #-----------------------------------------------------------------#
        Psi[0]  = np.gradient(Tau[2],h, edge_order=2)[1] -\
                    np.gradient(Tau[1], h, edge_order=2)[2]
        #-----------------------------------------------------------------#
        # Psi_{12}                                                        #
        #-----------------------------------------------------------------#
        Psi[1]  = np.gradient(Tau[4],h, edge_order=2)[1] -\
                    np.gradient(Tau[3], h, edge_order=2)[2]
        #-----------------------------------------------------------------#
        # Psi_{13}                                                        #
        #-----------------------------------------------------------------#
        Psi[2]  = np.gradient(Tau[5],h, edge_order=2)[1] -\
                    np.gradient(Tau[4], h, edge_order=2)[2]
        #-----------------------------------------------------------------#
        # Psi_{21}                                                        #
        #-----------------------------------------------------------------#
        Psi[3]  = np.gradient(Tau[0],h, edge_order=2)[2] -\
                    np.gradient(Tau[2], h, edge_order=2)[0]
        #-----------------------------------------------------------------#
        # Psi_{22}                                                        #
        #-----------------------------------------------------------------#
        Psi[4]  = np.gradient(Tau[1],h, edge_order=2)[2] -\
                    np.gradient(Tau[4], h, edge_order=2)[0]
        #-----------------------------------------------------------------#
        # Psi_{23}                                                        #
        #-----------------------------------------------------------------#
        Psi[5]  = np.gradient(Tau[2],h, edge_order=2)[2] -\
                    np.gradient(Tau[5], h, edge_order=2)[0]
        #-----------------------------------------------------------------#
        # Psi_{31}                                                        #
        #-----------------------------------------------------------------#
        Psi[6]  = np.gradient(Tau[1],h, edge_order=2)[0] -\
                    np.gradient(Tau[0], h, edge_order=2)[1]
        #-----------------------------------------------------------------#
        # Psi_{32}                                                        #
        #-----------------------------------------------------------------#
        Psi[7]  = np.gradient(Tau[3],h, edge_order=2)[0] -\
                    np.gradient(Tau[1], h, edge_order=2)[1]
        #-----------------------------------------------------------------#
        # Psi_{33}                                                        #
        #-----------------------------------------------------------------#
        Psi[8]  = np.gradient(Tau[4],h, edge_order=2)[0] -\
                    np.gradient(Tau[2], h, edge_order=2)[1]

    return Psi
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

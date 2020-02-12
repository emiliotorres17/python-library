#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the SGS production in
    the enstrophy transport equation.

    **** Note: The spectral differentiation still needs to be added

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
def SGS_production(
        w1,                 # vorticity-1 component
        w2,                 # vorticity-2 component
        w3,                 # vorticity-3 component
        Tau,                # SGS; (6,64,64,64)
        h       = False,    # spatial step size
        flag    = True):    # spectral flag; default gradient tool

    """ Calculating the SGS enstrophy production """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if h is False:
        Pi  = np.pi
        num = 64
        h   = (2.0*Pi)/num
    #---------------------------------------------------------------------#
    # Calculating psi                                                     #
    #---------------------------------------------------------------------#
    Psi = psi(Tau, h, flag)
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    dim         = np.shape(w1)[0]
    SGS_prod    = np.zeros((dim, dim, dim))
    #---------------------------------------------------------------------#
    # Calculating gradients                                               #
    #---------------------------------------------------------------------#
    grad_w1     = np.gradient(w1, h, edge_order=2)
    grad_w2     = np.gradient(w2, h, edge_order=2)
    grad_w3     = np.gradient(w3, h, edge_order=2)
    #---------------------------------------------------------------------#
    # Calculating SGS production                                          #
    #---------------------------------------------------------------------#
    if flag is False:               # spectral method
        print("Spectral differentiation has not been added")
        sys.exit(1)
    else:                           # gradient tool method
        #-----------------------------------------------------------------#
        # terms 1-3 (i=1)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[0]*grad_w1[0]
        SGS_prod    += Psi[1]*grad_w1[0]
        SGS_prod    += Psi[2]*grad_w1[0]
        #-----------------------------------------------------------------#
        # terms 4-6 (i=2)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[3]*grad_w2[0]
        SGS_prod    += Psi[4]*grad_w2[0]
        SGS_prod    += Psi[5]*grad_w2[0]
        #-----------------------------------------------------------------#
        # terms 7-9 (i=3)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[6]*grad_w3[0]
        SGS_prod    += Psi[7]*grad_w3[0]
        SGS_prod    += Psi[8]*grad_w3[0]

    return SGS_prod
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

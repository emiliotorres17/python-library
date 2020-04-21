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
def term(
        nu_sgs,         # strain rates
        b,              # B term in the enstrophy transport equation
        d):             # D term in the enstrophy transport equation

    """ Calculating term 1 in the SGS enstrophy study """
    #---------------------------------------------------------------------#
    # Term 1                                                              #
    #---------------------------------------------------------------------#
    term    = np.zeros((64,64,64))
    term    += (nu_sgs/nu)*(b+d)

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

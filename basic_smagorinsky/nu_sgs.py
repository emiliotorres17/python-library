#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate the turbulent viscosity for
    a Basic Smagorinsky simulation.

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
def nu_sgs_BS(
        Sij,
        cs      = True,         # C_{S} value; default C_{S} = 0.17
        delta   = True):        # Delta value; default is 2*pi/30

    """ Calculating the turbulent viscosity """
    #---------------------------------------------------------------------#
    # Setting default values                                              #
    #---------------------------------------------------------------------#
    if cs is True:
        cs  = 0.17
    if delta is True:
        delta = (2.0*np.pi)/30.0
    #---------------------------------------------------------------------#
    # Finding the magnitude of S_{ij}                                     #
    #---------------------------------------------------------------------#
    Smag    = np.zeros((64,64,64))
    Smag    += np.sqrt(2.0*(Sij[0]*Sij[0]+2.0*Sij[1]*Sij[1]+\
                        2.0*Sij[2]*Sij[2]+ Sij[3]*Sij[3]+2.0*Sij[4]*Sij[4]\
                        +Sij[5]*Sij[5]))
    #---------------------------------------------------------------------#
    # Calculating nu_{sgs}                                                #
    #--------------------------------------------------------------------#
    nu_sgs = -2.0*(cs*delta)**2.0*Smag

    return nu_sgs
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])

    print('**** This has not been unit tested ****')

    sys.exit(0)

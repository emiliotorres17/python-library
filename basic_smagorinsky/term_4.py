#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate term-4 in the enstrophy
    study equation, namely,

    epsilon_{ijk}*omega_{i}*S_{kl}*d2(nu_{sgs})/dx_{j}d_{j}

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
# term 4                                                                  #
#-------------------------------------------------------------------------#
def term_4(
        omega1,
        omega2,
        omega3,
        S,
        nu_sgs):

    """ Calculating term 4 """
    #---------------------------------------------------------------------#
    # defining domain variables and preallocating variables               #
    #---------------------------------------------------------------------#
    term        = np.zeros((64,64,64))
    h           = 2.0*np.pi/64.0
    #---------------------------------------------------------------------#
    # i=1, j=2, and k=3                                                   #
    #---------------------------------------------------------------------#
    term    += omega1*(S[2]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[1] +\
                        S[4]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[1] +\
                        S[5]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[1])
    #---------------------------------------------------------------------#
    # i=1, j=3, and k=2                                                   #
    #---------------------------------------------------------------------#
    term    -= omega1*(S[1]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[0] +\
                        S[3]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[0] +\
                        S[4]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[0])
    #---------------------------------------------------------------------#
    # i=2, j=3, and k=1                                                   #
    #---------------------------------------------------------------------#
    term    += omega2*(S[0]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[0] +\
                        S[1]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[0] +\
                        S[2]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[0])
    #---------------------------------------------------------------------#
    # i=2, j=1, and k=3                                                   #
    #---------------------------------------------------------------------#
    term    -= omega2*(S[2]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[1] +\
                        S[4]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[1] +\
                        S[5]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[1])
    #---------------------------------------------------------------------#
    # i=3, j=1, and k=2                                                   #
    #---------------------------------------------------------------------#
    term    += omega3*(S[1]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[1] +\
                        S[3]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[1] +\
                        S[4]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[1])
    #---------------------------------------------------------------------#
    # i=3, j=2, and k=1                                                   #
    #---------------------------------------------------------------------#
    term    -= omega3*(S[0]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[2], h, edge_order=2)[1] +\
                        S[1]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[1], h, edge_order=2)[1] +\
                        S[2]*np.gradient(\
                            np.gradient(nu_sgs, h, edge_order=2)[0], h, edge_order=2)[1])

    return term
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])

    print('**** This has not been unit tested ****')

    sys.exit(0)

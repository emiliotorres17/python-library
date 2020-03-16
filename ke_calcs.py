#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of these subroutines is to calculate the following:
        1. average kinetic energy
        2. kinetic energy field

    **** Note: This subroutine has not been unit tested, however a visual
            test has been performed between python and MATLAB in:
                ~/Projects/Fall-2019/MAE-792/blowup/test-10-2019/data-65-8

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
# Kinetic energy                                                          #
#-------------------------------------------------------------------------#
def ke_field(
        u1,                     # velocity-1 component
        u2,                     # velocity-2 component
        u3):                    # velocity-3 component

    """ Calculating the kinetic energy field """
    #---------------------------------------------------------------------#
    # Determining the shape of the KE field                               #
    #---------------------------------------------------------------------#
    dim     = u1.shape
    ke      = np.empty((dim[0], dim[1], dim[2], dim[3]))
    #---------------------------------------------------------------------#
    # Time loop                                                           #
    #---------------------------------------------------------------------#
    print_count = 0
    for i in range(0, dim[3]):
        term1           = np.square(u1[:,:,:,i])
        term2           = np.square(u2[:,:,:,i])
        term3           = np.square(u3[:,:,:,i])
        ke[:,:,:,i]     = 0.5*(term1 + term2 + term3)
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('KE field ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return ke
#-------------------------------------------------------------------------#
# Average kinetic energy                                                  #
#-------------------------------------------------------------------------#
def ke_average(
        u1,                 # velocity-1 component
        u2,                 # velocity-2 component
        u3):                # velocity-3 component

    """ Calculating the average kinetic energy """
    #---------------------------------------------------------------------#
    # Calculated average kinetic energy                                   #
    #---------------------------------------------------------------------#
    dim         = u1.shape
    keavg       = np.zeros(dim[-1])
    print_count = 0
    for i in range(0, dim[-1]):
        term1           = np.square(u1[:,:,:,i])
        term2           = np.square(u2[:,:,:,i])
        term3           = np.square(u3[:,:,:,i])
        ke              = 0.5*(term1 + term2 + term3)
        keavg[i]        = np.mean(ke)
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('KE average ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return keavg
#-------------------------------------------------------------------------#
# Average kinetic energy                                                  #
#-------------------------------------------------------------------------#
def ke_average2(
        KE):                # kinetic energy field

    """ Calculating the average kinetic energy """
    #---------------------------------------------------------------------#
    # Calculated average kinetic energy                                   #
    #---------------------------------------------------------------------#
    dim         = KE.shape
    keavg       = np.zeros(dim[-1])
    print_count = 0
    for i in range(0, dim[-1]):
        keavg[i]        = np.mean(KE[:,:,:,i])
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('KE average ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return keavg
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Print statement                                                     #
    #---------------------------------------------------------------------#
    print('This has not been unit tested')

    sys.exit(0)

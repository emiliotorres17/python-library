#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to locate the position of a particle
    in periodic cube.

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
# User defined function                                                   #
#=========================================================================#
#-------------------------------------------------------------------------#
# Cube locator                                                            #
#-------------------------------------------------------------------------#
def cube_locator(
        x):                     # position value

    """ Calculating the location/position in a periodical cube """
    #---------------------------------------------------------------------#
    # Defining needed variables                                           #
    #---------------------------------------------------------------------#
    pi      = np.pi
    #=====================================================================#
    # Periodic boundary condition locators                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Locations greater then 2.0*pi                                       #
    #---------------------------------------------------------------------#
    if x > 2.0*pi:
        N       = x/(2.0*pi)
        part    = N - np.floor(N)
        pos     = 2.0*pi*part
    #---------------------------------------------------------------------#
    # Locations from 2.0*pi <= x <= 0.0                                   #
    #---------------------------------------------------------------------#
    elif x < 0.0 and abs(x) < 2.0*pi:
        pos     = 2.0*pi - abs(x)
    #---------------------------------------------------------------------#
    # Locations from x < 2.0*pi                                           #
    #---------------------------------------------------------------------#
    elif x < 0.0 and abs(x)> 2.0*pi:
        N       = abs(x/(2.0*pi))
        part    = N - np.floor(N)
        pos     = 2.0*pi - 2.0*pi*part
    elif x > 0.0 and x < 2.0*pi:
        pos     = x
    else:
        print('**** Cannot find cube location x = %.5f ****'        %(x))
        sys.exit(1)

    return pos
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Defining testing variables                                          #
    #---------------------------------------------------------------------#
    xpos        = np.array([6.30, 500.25, 185.03, 7.20, 6.55])
                                        # test values x > 2.0*pi
    xcube       = np.array([5.83, 0.50, 4.13, 6.20, 5.11])
                                        # test values 0.0 <= x <= 2.0*pi
    xneg        = np.array([-0.23, -1.58, -3.65, -0.85, -7.85, -500.29])
                                        # test values x < 0.0
    solpos      = np.array([0.01681469282, 3.878360733, 2.817626092,\
                                0.9168146928, 0.2668146928])
                                        # solutions x > 2.0*pi
    solcube     = np.array([5.83, 0.50, 4.13, 6.20, 5.11])
                                        # solutions 0.0 <= x <= 2.0*pi
    solneg      = np.array([6.053185307, 4.703185307, 2.633185307,\
                                5.433185307, 4.716370614, 2.364824574])
                                        # solution for negative values
    tol         = 1e-05
    #---------------------------------------------------------------------#
    # Testing the positive values greater then 2.0*pi                     #
    #---------------------------------------------------------------------#
    print('**** Testing values greater then 2.0*pi ****')
    for i, sol in enumerate(solpos, 0):
        calc        = cube_locator(xpos[i])
        err         = abs((calc - sol)/sol)
        if err < tol:
            results = "***** Test #%i passed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)
        else:
            results = "***** Test #%i failed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)
    #---------------------------------------------------------------------#
    # Testing values between 0.0 <= x <= 2.0*pi                           #
    #---------------------------------------------------------------------#
    print('\n\n**** Testing values between 0.0 <= x <= 2.0*pi ****')
    for i, sol in enumerate(solcube, 0):
        calc        = cube_locator(xcube[i])
        err         = abs((calc - sol)/sol)
        if err < tol:
            results = "***** Test #%i passed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)
        else:
            results = "***** Test #%i failed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)
    #---------------------------------------------------------------------#
    # Testing values < 0.0                                                #
    #---------------------------------------------------------------------#
    print('\n\n**** Testing values x < 0.0 ****')
    for i, sol in enumerate(solneg, 0):
        calc        = cube_locator(xneg[i])
        err         = abs((calc - sol)/sol)
        if err < tol:
            results = "***** Test #%i passed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)
        else:
            results = "***** Test #%i failed \t\t error= = %.5e \t\t calc = %.7f \t\t exp = %.7f"\
                                    %(i + 1, err, calc, sol)
            print(results)

    sys.exit(0)

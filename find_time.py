#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to function that can find a specific
    value in a 1D vector.

    **** Note: This is an inclusive counter if the value is not present the
    subroutine returns the location of the first value passed the requested
    value.

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
# Finding a specific time value                                           #
#-------------------------------------------------------------------------#
def find_time(
        val,                    # value of desire
        Times):                 # time vector Nx1

    """ Finding a given time value for an time advancing vector """
    #---------------------------------------------------------------------#
    # Finding a specific time value                                       #
    #---------------------------------------------------------------------#
    for count, t in enumerate(Times):
        if t > val:
            Index = count - 1
            break
        #-----------------------------------------------------------------#
        # If the value is > then requested value the last value is given  #
        #-----------------------------------------------------------------#
        if count == len(Times)-1:
            Index = count

    return Index
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
    time    = np.linspace(0.0, 100, 501)
    time    = np.round(time, 3)
    #---------------------------------------------------------------------#
    # Test#1                                                              #
    #---------------------------------------------------------------------#
    index   = find_time(0.6, time)
    if index == 3:
        print("**** Test#1 passed \t\t expected = 3 \t\t calc = %i" \
                            %(index))
    else:
        print("**** Test#1 failed \t\t expected = 3 \t\t calc = %i"\
                            %(index))
        sys.exit(1)
    #---------------------------------------------------------------------#
    # Test#2                                                              #
    #---------------------------------------------------------------------#
    index   = find_time(0.8, time)
    if index == 4:
        print("**** Test#2 passed \t\t expected = 4 \t\t calc = %i" \
                            %(index))
    else:
        print("**** Test#2 failed \t\t expected = 4 \t\t calc = %i"\
                            %(index))
        sys.exit(1)
    #---------------------------------------------------------------------#
    # Test#3                                                              #
    #---------------------------------------------------------------------#
    index   = find_time(200.0, time)
    if index == 500:
        print("**** Test#3 passed \t\t expected = 500 \t\t calc = %i" \
                            %(index))
    else:
        print("**** Test#2 failed \t\t expected = 500 \t\t calc = %i"\
                            %(index))
        sys.exit(1)

    print('\n**** Successful Run ****')
    sys.exit(0)

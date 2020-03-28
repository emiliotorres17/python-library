#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of these subroutines are to find the minimum and maximum
    values and locations of an 2,3,4D numpy array.

    **** Still need to test the 4D arrays ****

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
# Finding the maximum subroutine 2D                                       #
#-------------------------------------------------------------------------#
def find_max2D(
        field):                 # 2D field

    """ Finding the maximum of 2D field """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmax(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Val     = np.amax(field)
    if Val != field[Loc1, Loc2]:
        print("**** Error in finding the maximum value")
        sys.exit(1)

    return Val, Loc1, Loc2
#-------------------------------------------------------------------------#
# Finding the maximum subroutine 3D                                       #
#-------------------------------------------------------------------------#
def find_max3D(
        field):                 # 3D field

    """ Finding the maximum of 3D field """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmax(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Loc3    = index[2]
    Val     = np.amax(field)
    if Val != field[Loc1, Loc2, Loc3]:
        print("**** Error in finding the maximum value")
        sys.exit(1)

    return Val, Loc1, Loc2, Loc3
#-------------------------------------------------------------------------#
# Finding the minimum subroutine 2D                                       #
#-------------------------------------------------------------------------#
def find_min2D(
        field):                 # 2D field

    """ Finding the minimum of 2D field """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmin(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Val     = np.amin(field)
    if Val != field[Loc1, Loc2]:
        print("**** Error in finding the minimum value")
        sys.exit(1)

    return Val, Loc1, Loc2
#-------------------------------------------------------------------------#
# Finding the minimum subroutine 3D                                       #
#-------------------------------------------------------------------------#
def find_min3D(
        field):                 # 3D field

    """ Finding the minimum of 3D field """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmin(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Loc3    = index[2]
    Val     = np.amin(field)
    if Val != field[Loc1, Loc2, Loc3]:
        print("**** Error in finding the minimum value")
        sys.exit(1)

    return Val, Loc1, Loc2, Loc3
#-------------------------------------------------------------------------#
# Finding the maximum subroutine 4D                                       #
#-------------------------------------------------------------------------#
def find_max4D(
        field):                 # 4D field

    """ Finding the maximum of 4D field """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmax(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Loc3    = index[2]
    Loc4    = index[3]
    Val     = np.amax(field)
    if Val != field[Loc1, Loc2, Loc3, Loc4]:
        print("**** Error in finding the maximum value")
        sys.exit(1)

    return Val, Loc1, Loc2, Loc3, Loc4
#-------------------------------------------------------------------------#
# Finding the minimum subroutine 4D                                       #
#-------------------------------------------------------------------------#
def find_min4D(
        field):                 # 4D field

    """ Finding the minimum of 4D field """
    #---------------------------------------------------------------------#
    # Finding the minimum subroutine                                      #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmin(field, axis=None), field.shape)
    Loc1    = index[0]
    Loc2    = index[1]
    Loc3    = index[2]
    Loc4    = index[3]
    Val     = np.amin(field)
    if Val != field[Loc1, Loc2, Loc3, Loc4]:
        print("**** Error in finding the minimum value")
        sys.exit(1)

    return Val, Loc1, Loc2, Loc3, Loc4
#-------------------------------------------------------------------------#
# Finding the maximum subroutine 4D                                       #
#-------------------------------------------------------------------------#
def find_max_per_time4D(
        field,                      # field variable
        name = 'kinetic energy'):   # name of the field

    """ Subroutine to find the maximum value per time step """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine 4D                                   #
    #---------------------------------------------------------------------#
    N           = field.shape[-1]
    max_vec     = np.zeros(N) 
    print_count = 0
    for i in range(0,N):
        max_vec[i] = np.amax(field[:,:,:,i])
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print( name + '---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return max_vec
#-------------------------------------------------------------------------#
# Finding the maximum subroutine 4D                                       #
#-------------------------------------------------------------------------#
def find_min_per_time4D(
        field,                      # field variable
        name = 'kinetic energy'):   # name of the field

    """ Subroutine to find the maximum value per time step """
    #---------------------------------------------------------------------#
    # Finding the maximum subroutine 4D                                   #
    #---------------------------------------------------------------------#
    N           = field.shape[-1]
    min_vec     = np.zeros(N) 
    print_count = 0
    for i in range(0,N):
        min_vec[i] = np.amin(field[:,:,:,i])
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print( name + '---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return min_vec
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Testing the maximum function 2D                                     #
    #---------------------------------------------------------------------#
    test                = np.array([[1,2,3], [-15,15,8]])
    sol                 = 15
    (val, loc1, loc2)   = find_max2D(test)
    if val == sol:
        print("**** Maximum test passed 2D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
        print("\tLocation 1: expected = 1 \t\t calculated = %i"\
                        %(loc1))
        print("\tLocation 2: expected = 1 \t\t calculated = %i"\
                        %(loc2))
    else:
        print("**** Maximum test failed 2D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
    #---------------------------------------------------------------------#
    # Testing the minimum function 3D                                     #
    #---------------------------------------------------------------------#
    test                = np.array([[1,2,3], [-15,15,8]])
    sol                 = -15
    (val, loc1, loc2)   = find_min2D(test)
    if val == sol:
        print("\n**** Minimum test passed 2D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
        print("\tLocation 1: expected = 1 \t\t calculated = %i"\
                        %(loc1))
        print("\tLocation 2: expected = 0 \t\t calculated = %i"\
                        %(loc2))
    else:
        print("**** Minimum test failed 2D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
    #---------------------------------------------------------------------#
    # Testing the maximum function 3D                                     #
    #---------------------------------------------------------------------#
    test                    = np.zeros((10,10,10))
    test[8,5,7]             = 8500
    sol                     = test[8,5,7]
    (val, loc1, loc2, loc3) = find_max3D(test)
    if val == sol:
        print("\n**** Maximum test passed 3D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
        print("\tLocation 1: expected = 8 \t\t calculated = %i"\
                        %(loc1))
        print("\tLocation 2: expected = 5 \t\t calculated = %i"\
                        %(loc2))
        print("\tLocation 3: expected = 7 \t\t calculated = %i"\
                        %(loc3))
    else:
        print("**** Maximum test failed 3D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
    #---------------------------------------------------------------------#
    # Testing the minimum function 3D                                     #
    #---------------------------------------------------------------------#
    test                    = np.zeros((10,10,10))
    test[6,4,1]             = -8500
    sol                     = test[6,4,1]
    (val, loc1, loc2, loc3) = find_min3D(test)
    if val == sol:
        print("\n**** Minimum test passed 3D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))
        print("\tLocation 1: expected = 6 \t\t calculated = %i"\
                        %(loc1))
        print("\tLocation 2: expected = 4 \t\t calculated = %i"\
                        %(loc2))
        print("\tLocation 3: expected = 1 \t\t calculated = %i"\
                        %(loc3))
    else:
        print("**** Minimum test failed 3D \t\t expected = %i \t\t calculated = %i"\
                        %(sol, val))

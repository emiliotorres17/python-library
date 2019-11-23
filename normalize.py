#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to perform various types of
    normalizations for both 2D and 3D fields.

    **** More extensive unit testing still needs to be done ****

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
from min_max    import find_max2D
from min_max    import find_min2D
from min_max    import find_max3D
from min_max    import find_min3D
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Normalizing from 0 to 1                                                 #
#-------------------------------------------------------------------------#
def norm_min_max(
        field):                 # field that is being normalized

    """ Normalizing a field by dividing subtracting the minimum and
    dividing by the maximum, so all the values are between 0 <= x <= 1 """
    #---------------------------------------------------------------------#
    # Finding the minimum and maximum                                     #
    #---------------------------------------------------------------------#
    dim     = field.shape
    if len(dim) == 2:
        minval  = find_min2D(field)[0]
        maxval  = find_max2D(field)[0]
    elif len(dim) == 3:
        minval  = find_min3D(field)[0]
        maxval  = find_max3D(field)[0]
    else:
        print('Dimensions out out of bound')
        print('\tDimensions =' + str(dim))
        sys.exit(1)
    #---------------------------------------------------------------------#
    # Normalizing                                                         #
    #---------------------------------------------------------------------#
    field       = field - minval
    field       = field/(maxval - minval)

    return field
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(["clear"])
    #=====================================================================#
    # Simple unit testing (Need more cases)                               #
    #=====================================================================#
    a           = np.array([[1,2,3], [2,3,5,], [4,3,5]])
    anorm       = norm_min_max(a)
    valmin      = find_min2D(anorm)[0]
    valmax      = find_max2D(anorm)[0]
    #---------------------------------------------------------------------#
    # Printing the normalize array for a visual check                     #
    #---------------------------------------------------------------------#
    print('Input array:')
    print(a)
    print('\nOutput array:')
    print(anorm)
    #---------------------------------------------------------------------#
    # Testing the minimum                                                 #
    #---------------------------------------------------------------------#
    if valmin == 0.0:
        print('\n\n**** Test passed minimum check')
        print('\t minimum --> expected = 0.0 \t calc = %.1f'\
                        %(valmin))
    else:
        print('**** Test failed minimum check')
        print('\t minimum --> expected = 0.0 \t calc = %.1f'\
                        %(valmin))
        sys.exit(1)
    #---------------------------------------------------------------------#
    # Testing the maximum                                                 #
    #---------------------------------------------------------------------#
    if valmax == 1.0:
        print('\n**** Test passed maximum check')
        print('\t maximum --> expected = 1.0 \t calc = %.1f'\
                        %(valmax))
    else:
        print('\n**** Test failed maximum check')
        print('\t maximum --> expected = 1.0 \t calc = %.1f'\
                        %(valmax))
        sys.exit(1)

    sys.exit(0)

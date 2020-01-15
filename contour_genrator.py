#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to plot any quantity from LES pseudo
    spectral code.

    **** Note: Makes subplots for a 64x64x64 quantity with the x_{1} and
               the x_{2} domain ranging from 0<= x_{i} < = 2.0*pi. This
               domain has been hard coded into the subroutine, however
               adjusting the code for a different domain should be pretty
               straight forward.

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
import numpy as np
import matplot.pyplot as plt
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Contour subplots                                                        #
#-------------------------------------------------------------------------#
def subplot_gen(
        var,                    # variable
        name,                   # file name
        loc,                    # location of where to save
        color   = "jet",        # color scheme (default set to jet)
        cb      = False,        # color bar flag (default set to off)
        xlab    = False,        # x-label flag (default set to off)
        ylab    = False):       # y-label flag (default set to off)

    """ Subroutine to generate subplots for any quantity (64x64x64) """
    #---------------------------------------------------------------------#
    # Checking the location path to make sure it has a OS separator       #
    #---------------------------------------------------------------------#
    sep     = os.sep                        # OS separator
    if loc[-1] != sep:                      # adding the separator
        loc     = loc + sep
    #---------------------------------------------------------------------#
    # Verifying the location                                              #
    #---------------------------------------------------------------------#
    isFile  = os.path.isdir(loc)            # checking existence of path
    if isFile is False:
        os.mkdir(loc)                       #  creating directory
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 64                                # number of data pts.
    pi          = np.pi                             # pi
    X1          = np.linspace(0.0, 2.0*pi, N)       # x_{1} domain
    X2          = np.linspace(0.0, 2.0*pi, N)       # x_{2} domain
    (X1, X2)    = np.meshgrid(X1, X2)               # mesh grid (plotting)
    #---------------------------------------------------------------------#
    # Generating contours                                                 #
    #---------------------------------------------------------------------#
    print_count     = 20                    # print counter
    dim             = var.shape             # dimensions of the quantity
    for i in range(0, dim[2]):
        plt.subplot(8,8,i+1)
        cnt     = plt.contourf(X1, X2, var[:,:,i], 500, cmap=color)
        for c in cnt.collections:           # setting the edge colors
            c.set_edgecolors("faces")
        if cb is not False:                 # generating the color bar
            plt.colorbar()
        if xlab is not False:               # generating the x_{1}-label
            plt.xlabel("$0 \leq x_{1} \leq 2\pi$")
        if ylab is not False:               # generating the x_{2}-label
            plt.ylabel("$0 \leq x_{1} \leq 2\pi$")
        #-----------------------------------------------------------------#
        # Print counter                                                   #
        #-----------------------------------------------------------------#
        if print_count > 10:
            print("z step --> %i"               %(i))
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Saving and clearing figure                                          #
    #---------------------------------------------------------------------#
    plt.savefig(loc + name + ".pdf")        # storing figure
    plt.clf()                               # clearing figure

    return
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Warning display                                                     #
    #---------------------------------------------------------------------#
    print("***** Warning this subroutine has NOT been tested ****")
    sys.exit(0)

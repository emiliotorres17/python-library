#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is perform a variety of different order
    interpolation schemes.

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
import matplotlib.pyplot as plt
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Find x values                                                           #
#-------------------------------------------------------------------------#
def find_xval(
        xvec,
        xval):

    """ returns grid coordinate of the n+1 of the value; xvec has to start
        at zero or higher """
    #---------------------------------------------------------------------#
    # Looping over x values                                               #
    #---------------------------------------------------------------------#
    for i in range(0, len(xvec)):
        if xvec[i] > xval:
            break
    print(xvec[i])
    return i
#-------------------------------------------------------------------------#
# 2nd order 1D interpolation                                              #
#-------------------------------------------------------------------------#
def interp_1D_2(
        x,                  # x values
        y,                  # y values
        xval):              # x pt. of interest

    """ Calculating a point based off of 2nd order linear interpolation """
    #---------------------------------------------------------------------#
    # Interpolating                                                       # 
    #---------------------------------------------------------------------#
    coor    = find_xval(x,xval)
    m       = (y[coor] - y[coor-1])/(x[coor] - x[coor-1])
    yval    = y[coor] - m*(x[coor]-xval)   

    return yval
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    N   = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, \
                int(4096*2)]
    err = np.zeros(len(N))
    dx  = np.zeros(len(N))
    #---------------------------------------------------------------------#
    # Looping over N                                                      #
    #---------------------------------------------------------------------#
    for i, n in enumerate(N):
        print(n)
        xint    = 0.554654512121354774524567578
        xvec    = np.linspace(0.0, 1.0, n+1)
        yvec    = np.sin(xvec)**2.0 + 3.0*np.cos(xvec)**2.0
        ycalc   = interp_1D_2(xvec, yvec, xint)  
        yexac   = np.sin(xint)**2.0 + 3.0*np.cos(xint)**2.0
        err[i]  = abs(yexac -  ycalc) 
        dx[i]   = 1.0/n
    #---------------------------------------------------------------------#
    # Font settings                                                       #
    #---------------------------------------------------------------------#
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    #SMALL_SIZE = 14
    #MEDIUM_SIZE = 18
    #BIGGER_SIZE = 12
    #plt.rc('font',      size=SMALL_SIZE)            # controls default text sizes
    #plt.rc('axes',      titlesize=SMALL_SIZE)       # fontsize of the axes title
    #plt.rc('axes',      labelsize=MEDIUM_SIZE)      # fontsize of the x and y labels
    #plt.rc('xtick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    #plt.rc('ytick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    #plt.rc('legend',    fontsize=SMALL_SIZE)        # legend fontsize
    #plt.rc('figure',    titlesize=BIGGER_SIZE)      # fontsize of the figure title
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    plt.loglog(dx, err, 'ro--', lw=1.5)
    plt.grid(True)
    #plt.yticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2])
    #plt.yticks(np.linspace(1e-16, 1e-2, 10))
    #plt.ylim([1e-16, 1e-02])
    plt.show()

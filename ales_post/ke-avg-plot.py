#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to plot the average kinetic energy for
    the different CDS values.

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Purpose                                                                 #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from ales_post.plot_settings        import plot_setting
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    media_path  = pwd + '%c..%cmedia%c'             %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    ts      = 0
    tf      = -1
    data    = loadtxt('../ke-omega-avg.txt', skiprows=2)
    steps   = data[ts:tf,0]
    time    = data[ts:tf,1]-15.0
    ke      = data[ts:tf,2]
    enst    = data[ts:tf,3]
    print(ke[-1])
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
    #---------------------------------------------------------------------#
    plot_setting()
    fig, [ax1, ax2]     = plt.subplots(2,1)
    #---------------------------------------------------------------------#
    # Plotting average kinetic energy v. time                             #
    #---------------------------------------------------------------------#
    ax1.plot([0,0], [amin(ke), amax(ke)], 'b--', lw=1.5)
    ax1.plot(time, ke, 'r', lw=1.5) 
    ax1.set_xlabel('Time')
    ax1.set_ylabel('$\\langle k \\rangle$')
    ax1.grid(True)
    #---------------------------------------------------------------------#
    # Plotting average enstrophy v. time                                  #
    #---------------------------------------------------------------------#
    ax1.plot([0,0], [amin(enst), amax(enst)], 'b--', lw=1.5)
    ax2.plot(time, enst, 'r', lw=1.5) 
    ax2.set_xlabel('Time')
    ax2.set_ylabel('$\\langle \\Omega \\rangle$')
    ax2.grid(True)
    plt.show()
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
    #---------------------------------------------------------------------#
    fig, [ax1, ax2]     = plt.subplots(2,1)
    #---------------------------------------------------------------------#
    # Plotting average kinetic energy v. time                             #
    #---------------------------------------------------------------------#
    ax1.plot(steps, ke, 'r', lw=1.5) 
    ax1.set_xlabel('Time steps')
    ax1.set_ylabel('$\\langle k \\rangle$')
    ax1.grid(True)
    #---------------------------------------------------------------------#
    # Plotting average enstrophy v. time                                  #
    #---------------------------------------------------------------------#
    ax2.plot(steps, enst, 'r', lw=1.5) 
    ax2.set_xlabel('Time steps')
    ax2.set_ylabel('$\\langle \\Omega \\rangle$')
    ax2.grid(True)
    plt.show()
    


    print('**** Successful run ****')
    sys.exit(0)

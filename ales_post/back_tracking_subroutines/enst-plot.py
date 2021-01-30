#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to plot the extracted values from the
    tracking subroutines.

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
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Plotting terms                                                          #
#-------------------------------------------------------------------------#
def plot_terms(
        DT,
        term,
        steps,
        sym,
        labs,
        name):

    """ Plotting the values of the different terms """
    #---------------------------------------------------------------------#
    # Plotting terms                                                      #
    #---------------------------------------------------------------------#
    plt.yscale('symlog')
    plt.plot([0, 0], [amin(DT), amax(DT)], 'b--', lw=1.5)
    plt.plot([-130, 200], [0, 0], 'b--', lw=1.5)
    if amax(DT) != 0:
        plt.plot(steps, DT, 'r', lw=1.5,\
                    label='$\\frac{1}{\\Omega}\\frac{D\\Omega}{Dt}$')
    else:
        plt.plot([0, 0], [-150, 150], 'b--', lw=1.5)
    plt.plot(steps, term, sym, lw=1.5, label=labs)
    plt.xlabel('Time step')
    plt.ylabel('$\\frac{1}{\\Omega} \\frac{D\\Omega}{Dt}$ Contributions')
    plt.grid(True)
    plt.legend(loc=2)
    plt.tight_layout()
    plt.xlim([-20, 20])
    plt.savefig(name, bbox_inches = "tight")
    #plt.show()
    plt.close()
    
    return
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
    media_path  = pwd + '%c..%cmedia%c'         %(sep, sep, sep)
    data_path   = pwd + '%c..%cdata%c'          %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Loading enstrophy data                                              #
    #---------------------------------------------------------------------#
    ts          = 0
    tint        = 240
    iters       = flipud(linspace(ts-222, tint-222, (tint-ts)+1))
    time        = loadtxt(data_path + 'enst.out', skiprows=1)[:,0]
    enst        = loadtxt(data_path + 'enst.out', skiprows=1)[:,-1]
    A_enst      = loadtxt(data_path + 'A-enst.out', skiprows=1)[:,-1]
    B_enst      = loadtxt(data_path + 'B-enst.out', skiprows=1)[:,-1]
    D_enst      = loadtxt(data_path + 'D-enst.out', skiprows=1)[:,-1]
    C_enst      = -loadtxt(data_path + 'C-enst.out', skiprows=1)[:,-1]
    P_enst      = -loadtxt(data_path + 'P-enst.out', skiprows=1)[:,-1]
    #---------------------------------------------------------------------#
    # Verification flag                                                   #
    #---------------------------------------------------------------------#
    ver_flag    = True
    if os.path.exists('enst-verification.out') is False:
        ver_flag    = False
    #---------------------------------------------------------------------#
    # Verification print                                                  #
    #---------------------------------------------------------------------#
    text        = ''
    for i in range(0, len(A_enst)):
        text    += '%10.1f\t%10.3f\n'               %(iters[i], time[i])            
    f   = open('enst-verification.out', 'w')
    f.write(text)
    f.close()
    if ver_flag is False:
        print('**** Check verification ****\n')
        sys.exit(8)
    #---------------------------------------------------------------------#
    # Material derivative of enstrophy                                    #
    #---------------------------------------------------------------------#
    dt_field    = (A_enst + B_enst + C_enst + D_enst + P_enst)/enst
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    print('Plotting:')
    plot_setting()
    plot_terms(dt_field, A_enst/enst, iters, 'b', '$A_{\\Omega}/\\Omega$',\
                    media_path + 'A-enst-240.png')
    print('\tA-enst')
    plot_terms(dt_field, C_enst/enst, iters, 'k', '$C_{\\Omega}/\\Omega$',\
                    media_path + 'C-enst-240.png')
    print('\tC-enst')
    plot_terms(dt_field, P_enst/enst, iters, 'g', '$P_{\\Omega}/\\Omega$',\
                    media_path + 'P-enst-240.png')
    print('\tP-enst')
    plot_terms(dt_field, (P_enst+A_enst+C_enst)/enst, iters, 'k',\
                    '$(P_{\\Omega}+C_{\\Omega} + A_{\\Omega})/\\Omega$',\
                    media_path + 'P-C-A-enst-240.png')
    print('\tP-C-A-enst')
    plot_terms(dt_field, (P_enst+C_enst)/enst, iters, 'k',\
                    '$(P_{\\Omega}+C_{\\Omega})/\\Omega$',\
                    media_path + 'P-C-enst-240.png')
    print('\tP-C-enst')
    plot_terms([0], (P_enst/C_enst), iters, 'k',\
                    '$P_{\\Omega}/C_{\\Omega}$',\
                    media_path + 'P-C-ratio-enst-240.png')
    print('\tP-C-ratio-enst')
    plot_terms(dt_field, B_enst/enst, iters, 'm', '$B_{\\Omega}/\\Omega$',\
                    media_path + 'B-enst-240.png')
    print('\tB-enst')
    plot_terms(dt_field, D_enst/enst, iters, 'c', '$D_{\\Omega}/\\Omega$',\
                    media_path + 'D-enst-240.png')
    print('\tD-enst')

    print('**** Successful run ****')
    sys.exit(0)

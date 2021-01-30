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
                    label='$\\frac{1}{k}\\frac{Dk}{Dt}$')
    else:
        plt.plot([0, 0], [-150, 150], 'b--', lw=1.5)
    plt.plot(steps, term, sym, lw=1.5, label=labs)
    plt.xlabel('Time step')
    plt.ylabel('$\\frac{1}{k} \\frac{Dk}{Dt}$ Contributions')
    plt.grid(True)
    plt.legend(loc=2)
    plt.tight_layout()
    plt.xlim([x_min_lim, x_max_lim])
    plt.savefig(name, bbox_inches = "tight")
    plt.show()
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
    # Loading ke data                                                     #
    #---------------------------------------------------------------------#
    ts          = 0
    tzero       = 215
    tint        = 255
    x_min_lim   = -int(tint-tzero)
    x_max_lim   = int(tint-tzero)
    iters       = flipud(linspace(ts-tzero, tint-tzero, (tint-ts)+1))
    time        = loadtxt(data_path + 'ke.out', skiprows=1)[:,0]
    ke          = loadtxt(data_path + 'ke.out', skiprows=1)[:,-1]
    A_ke        = loadtxt(data_path + 'A-ke.out', skiprows=1)[:,-1]
    B_ke        = loadtxt(data_path + 'B-ke.out', skiprows=1)[:,-1]
    D_ke        = loadtxt(data_path + 'D-ke.out', skiprows=1)[:,-1]
    C_ke        = loadtxt(data_path + 'C-ke.out', skiprows=1)[:,-1]
    P_ke        = loadtxt(data_path + 'P-ke.out', skiprows=1)[:,-1]
    #---------------------------------------------------------------------#
    # Verification flag                                                   #
    #---------------------------------------------------------------------#
    ver_flag    = True
    if os.path.exists('ke-verification.out') is False:
        ver_flag    = False
    #---------------------------------------------------------------------#
    # Verification print                                                  #
    #---------------------------------------------------------------------#
    text        = ''
    for i in range(0, len(A_ke)):
        text    += '%10.1f\t%10.3f\n'               %(iters[i], time[i])            
    f   = open('ke-verification.out', 'w')
    f.write(text)
    f.close()
    if ver_flag is False:
        print('**** Check verification ****\n')
        sys.exit(8)
    #---------------------------------------------------------------------#
    # Material derivative of enstrophy                                    #
    #---------------------------------------------------------------------#
    dt_field    = (A_ke + B_ke + C_ke + D_ke + P_ke)/ke
    #---------------------------------------------------------------------#
    # Plotting                                                            #
    #---------------------------------------------------------------------#
    print('Plotting:')
    plot_setting()
    plot_terms(dt_field, A_ke/ke, iters, 'b', '$A_{k}/k$',\
                    media_path + 'A-ke-240.png')
    print('\tA-ke')
    plot_terms(dt_field, C_ke/ke, iters, 'k', '$C_{k}/k$',\
                    media_path + 'C-ke-240.png')
    print('\tC-ke')
    plot_terms(dt_field, P_ke/ke, iters, 'g', '$P_{k}/k$',\
                    media_path + 'P-ke-25-240.png')

    print('\tP-ke')
    plot_terms(dt_field, (P_ke+A_ke+C_ke)/ke, iters, 'k',\
                    '$(P_{k}+C_{k} + A_{k})/k$',\
                    media_path + 'P-C-A-ke-240.png')
    print('\tP-C-A-ke')
    plot_terms(dt_field, (P_ke+C_ke)/ke, iters, 'k',\
                    '$(P_{k}+C_{k})/k$',\
                    media_path + 'P-C-ke-240.png')
    print('\tP-C-ke')
    plot_terms([0], (P_ke/C_ke), iters, 'k',\
                    '$P_{k}/C_{k}$',\
                    media_path + 'P-C-ratio-ke-240.png')
    print('\tP-C-ratio-ke')
    plot_terms(dt_field, B_ke/ke, iters, 'm', '$B_{k}/k$',\
                    media_path + 'B-ke-240.png')
    print('\tB-ke')
    plot_terms(dt_field, D_ke/ke, iters, 'c', '$D_{k}/k$',\
                    media_path + 'D-ke-240.png')
    print('\tD-ke')

    print('**** Successful run ****')
    sys.exit(0)

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
    data_path   = pwd + '%c..%cdata%c'          %(sep, sep, sep)
    media_path  = pwd + '%c..%cmedia%c'         %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    print('Loaded:')
    time        = load(data_path + 'time.npy')
    print('time')
    enst        = load(data_path + 'enst.npy')
    print('enst')
    ke          = load(data_path + 'ke.npy')
    print('ke')
    #---------------------------------------------------------------------#
    # Better output                                                       #
    #---------------------------------------------------------------------#
    text        = '%s'          %('step'.center(30))
    text        += '%s'         %('time'.center(30))
    text        += '%s'         %('<k>'.center(30))
    text        += '%s\n'       %('<Omega>'.center(30))
    count       = 0
    for i in range(0, len(time)): 
        count   += 1
        text    += '%s'                 %(str(i).center(30))
        text    += '%30.16e'            %(time[i])
        text    += '%30.16e'            %(mean(ke[:,:,:,i]))
        text    += '%30.16e\n'          %(mean(enst[:,:,:,i]))
        if count == 20:
            print('time step --> %i'        %(i))
            count   = 0 
    #---------------------------------------------------------------------#
    # Writting                                                            #
    #---------------------------------------------------------------------#
    f   = open('time-step-data.out', 'w')
    f.write(text)
    f.close()

    print('**** Successful run ****')
    sys.exit(0)

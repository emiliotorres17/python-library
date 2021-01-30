#!/usr/bin/env python3
"""========================================================================
Purpose:
    Back track the particle to compare to the forward tracking.

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
from numpy import *
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from ales_post.plot_settings    import plot_setting
from min_max                    import find_max3D
from rk4_back_tracking          import rk4_back_tracking
from rk4_back_tracking          import ghost_cells
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
    #---------------------------------------------------------------------#
    # Deleting verification files                                         #
    #---------------------------------------------------------------------#
    enst_name   = 'enst-verification.out'
    ke_name     = 'ke-verification.out'
    if os.path.exists(enst_name) is True:
        os.remove(enst_name)
    if os.path.exists(ke_name) is True:
        os.remove(ke_name)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    N       = 64
    dx      = 2.*pi/N
    x       = linspace(0.5*dx, 2.*pi-0.5*dx, N)
    y       = linspace(0.5*dx, 2.*pi-0.5*dx, N)
    z       = linspace(0.5*dx, 2.*pi-0.5*dx, N)
    #---------------------------------------------------------------------#
    # Loading enstrophy                                                   #
    #---------------------------------------------------------------------#
    print('Loading data:')
    ts          = 0 
    tf          = 251 
    tint        = 240
    time        = load(data_path + 'time.npy')[ts:tf]
    dt          = zeros(len(time)-1)
    count       = zeros(len(time)-1)
    print('\ttime')
    enst        = load(data_path + 'enst.npy')[:,:,:,ts:tf]    
    print(enst.shape)
    print('\tenstrophy')

    print(time[0])
    print(mean(enst[:,:,:,0]))
    
    print(time[-1])
    print(mean(enst[:,:,:,-1]))
    #---------------------------------------------------------------------#
    # Maximum values and coordinates                                      #
    #---------------------------------------------------------------------#
    [val, xc, yc, zc]   = find_max3D(enst[:,:,:,tint])
    print('enst max --> %10.5e'                 %(val))
    print('x-index = %3i\tx-value = %8.4e'      %(xc, x[xc]))
    print('y-index = %3i\ty-value = %8.4e'      %(yc, y[yc]))
    print('z-index = %3i\tz-value = %8.4e'      %(zc, z[zc]))
    #---------------------------------------------------------------------#
    # Initial point                                                       #
    #---------------------------------------------------------------------#
    x0  = 0.6652053603967715
    y0  = 5.0357656588257953
    z0  = 4.6776113491502915
    t0  = 15.3050502065754515
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    dx          = 2.*pi/N
    x1          = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    x2          = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    x3          = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    #---------------------------------------------------------------------#
    # Loading velocity                                                    #
    #---------------------------------------------------------------------#
    print('Loading data:')
    u1          = load(data_path + 'velocity1.npy')[:,:,:,ts:tf]
    u1          = ghost_cells(N, u1)
    print('\tvelocity-1')
    u2          = load(data_path + 'velocity2.npy')[:,:,:,ts:tf]
    u2          = ghost_cells(N, u2)
    print('\tvelocity-2')
    u3          = load(data_path + 'velocity3.npy')[:,:,:,ts:tf]
    u3          = ghost_cells(N, u3)
    print('\tvelocity-3')
    #---------------------------------------------------------------------#
    # Initial start point                                                 #
    #---------------------------------------------------------------------#
    print('%25.20f'     %(t0))
    print('%25.20f'     %(x0))
    print('%25.20f'     %(y0))
    print('%25.20f'     %(z0))
    #---------------------------------------------------------------------#
    # Back tracking                                                       #
    #---------------------------------------------------------------------#
    rk4_back_tracking([x0, y0, z0, t0], u1, u2, u3, \
                        x1, x2, x3, time, tint, 0,\
                        'enstrophy-back-print.out', 'enstrophy-back-coordinates.out')
    call(['python3', 'plot-coordinates.py'])
    print('**** Successful run ****')
    sys.exit(0)

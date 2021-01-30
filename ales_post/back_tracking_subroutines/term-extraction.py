#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to extract the transport terms using the
    coordinate values.

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
from scipy.interpolate  import RegularGridInterpolator as Reginterp
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from ales_post.plot_settings        import plot_setting
from rk4_back_tracking              import ghost_cells
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Extracting values of interest                                           #
#-------------------------------------------------------------------------#
def extract_values(
        func,
        X1,
        X2,
        X3,
        T, 
        name):

    """ subroutine to extract values """
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    val     = zeros(len(X1))
    output  = '%s %s %s %s %s\n'    %('time'.center(20), 'x1'.center(20),\
                                        'x2'.center(20), 'x3'.center(20),\
                                        'value'.center(55))
    #---------------------------------------------------------------------#
    # Extracting values of interest                                       #
    #---------------------------------------------------------------------#
    for i in range(0, len(X1)):
        print(i)
        val[i]  = func((X1[i], X2[i], X3[i], T[i])) 
        output  += '%20.16f %20.16f %20.16f %20.16f %25.16e\n'\
                        %(T[i], X1[i], X2[i], X3[i], val[i])
    #---------------------------------------------------------------------#
    # Writing data                                                        #
    #---------------------------------------------------------------------#
    f   = open(name, 'w')
    f.write(output)
    f.close()

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
    data_path   = pwd + '%c..%cdata%c'          %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    ts      = 0 
    tf      = 251
    N       = 64
    dx      = 2.*pi/N
    x       = linspace(-0.5*dx, 2.*pi+0.5*dx, N+2)
    y       = linspace(-0.5*dx, 2.*pi+0.5*dx, N+2)
    z       = linspace(-0.5*dx, 2.*pi+0.5*dx, N+2)
    #---------------------------------------------------------------------#
    # Extracting coordinate points                                        #
    #---------------------------------------------------------------------#
    coor    = loadtxt('enstrophy-back-coordinates.out', skiprows=1)
    t       = coor[:,0] 
    x1      = coor[:,1]
    x2      = coor[:,2]
    x3      = coor[:,3]
    #---------------------------------------------------------------------#
    # Flags                                                               #
    #---------------------------------------------------------------------#
    ke_enst_flag        = True
    ke_trans_flag       = True
    enst_trans_flag     = True
    cs2_flag            = False
    #---------------------------------------------------------------------#
    # KE                                                                  #
    #---------------------------------------------------------------------#
    if ke_trans_flag is True:
        #-----------------------------------------------------------------#
        # Loading data                                                    #
        #-----------------------------------------------------------------#
        print('Load data:')
        time    = load(data_path + 'time.npy')[ts:tf]
        print('\ttime')
        A       = load(data_path + 'A-ke.npy')[:,:,:,ts:tf]
        A       = ghost_cells(N, A) 
        print('\tA-ke')
        B       = load(data_path + 'B-ke.npy')[:,:,:,ts:tf]
        B       = ghost_cells(N, B) 
        print('\tB-ke')
        D       = load(data_path + 'D-ke.npy')[:,:,:,ts:tf]
        D       = ghost_cells(N, D) 
        print('\tD-ke')
        C       = load(data_path + 'C-ke.npy')[:,:,:,ts:tf]
        C       = ghost_cells(N, C) 
        print('\tC-ke')
        P       = load(data_path + 'P-ke.npy')[:,:,:,ts:tf]
        P       = ghost_cells(N, P)
        print('\tP-ke')
        #-----------------------------------------------------------------#
        # Generating interpolation functions                              #
        #-----------------------------------------------------------------#
        Afunc   = Reginterp((x,y,z,time,), A, method='linear') 
        Bfunc   = Reginterp((x,y,z,time,), B, method='linear') 
        Dfunc   = Reginterp((x,y,z,time,), D, method='linear') 
        Cfunc   = Reginterp((x,y,z,time,), C, method='linear') 
        Pfunc   = Reginterp((x,y,z,time,), P, method='linear') 
        #-----------------------------------------------------------------#
        # Extracting the values of interest                               #
        #-----------------------------------------------------------------#
        print('Extracting values of interest:')
        extract_values(Afunc, x1, x2, x3, t, data_path + 'A-ke.out')
        print('\tA-ke')
        extract_values(Bfunc, x1, x2, x3, t, data_path + 'B-ke.out')
        print('\tB-ke')
        extract_values(Dfunc, x1, x2, x3, t, data_path + 'D-ke.out')
        print('\tD-ke')
        extract_values(Cfunc, x1, x2, x3, t, data_path + 'C-ke.out')
        print('\tC-ke')
        extract_values(Pfunc, x1, x2, x3, t, data_path + 'P-ke.out')
        print('\tP-ke')
        #-----------------------------------------------------------------#
        # Deleting variables                                              #
        #-----------------------------------------------------------------#
        del A
        del B
        del D
        del C
        del P
        del Afunc
        del Bfunc
        del Dfunc
        del Cfunc
        del Pfunc
    #---------------------------------------------------------------------#
    # Enstrophy                                                           #
    #---------------------------------------------------------------------#
    if enst_trans_flag is True:
        #-----------------------------------------------------------------#
        # Loading data                                                    #
        #-----------------------------------------------------------------#
        print('Load data:')
        time    = load(data_path + 'time.npy')[ts:tf]
        print('\ttime')
        A       = load(data_path + 'A-enst.npy')[:,:,:,ts:tf]
        A       = ghost_cells(N, A)
        print('\tA-enst')
        B       = load(data_path + 'B-enst.npy')[:,:,:,ts:tf]
        B       = ghost_cells(N, B)
        print('\tB-enst')
        D       = load(data_path + 'D-enst.npy')[:,:,:,ts:tf]
        D       = ghost_cells(N, D)
        print('\tD-enst')
        C       = load(data_path + 'C-enst.npy')[:,:,:,ts:tf]
        C       = ghost_cells(N, C)
        print('\tC-enst')
        P       = load(data_path + 'P-enst.npy')[:,:,:,ts:tf]
        P       = ghost_cells(N, P)
        print('\tP-enst')
        #-----------------------------------------------------------------#
        # Generating interpolation functions                              #
        #-----------------------------------------------------------------#
        Afunc   = Reginterp((x,y,z,time,), A, method='linear') 
        Bfunc   = Reginterp((x,y,z,time,), B, method='linear') 
        Dfunc   = Reginterp((x,y,z,time,), D, method='linear') 
        Cfunc   = Reginterp((x,y,z,time,), C, method='linear') 
        Pfunc   = Reginterp((x,y,z,time,), P, method='linear') 
        #-----------------------------------------------------------------#
        # Extracting the values of interest                               #
        #-----------------------------------------------------------------#
        print('Extracting values of interest:')
        extract_values(Afunc, x1, x2, x3, t, data_path + 'A-enst.out')
        print('\tA-enst')
        extract_values(Bfunc, x1, x2, x3, t, data_path + 'B-enst.out')
        print('\tB-enst')
        extract_values(Dfunc, x1, x2, x3, t, data_path + 'D-enst.out')
        print('\tD-enst')
        extract_values(Cfunc, x1, x2, x3, t, data_path + 'C-enst.out')
        print('\tC-enst')
        extract_values(Pfunc, x1, x2, x3, t, data_path + 'P-enst.out')
        print('\tP-enst')
        #-----------------------------------------------------------------#
        # Deleting variables                                              #
        #-----------------------------------------------------------------#
        del A
        del B
        del D
        del C
        del P
        del Afunc
        del Bfunc
        del Dfunc
        del Cfunc
        del Pfunc
    #---------------------------------------------------------------------#
    # Enstrophy                                                           #
    #---------------------------------------------------------------------#
    if ke_enst_flag is True:
        #-----------------------------------------------------------------#
        # Loading data                                                    #
        #-----------------------------------------------------------------#
        print('Load data:')
        time    = load(data_path + 'time.npy')[ts:tf]
        print('\ttime')
        ke      = load(data_path + 'ke.npy')[:,:,:,ts:tf]
        ke      = ghost_cells(N, ke)
        print('\tke-field')
        enst    = load(data_path + 'enst.npy')[:,:,:,ts:tf]
        enst    = ghost_cells(N, enst)
        print('\tenst-enst')
        #-----------------------------------------------------------------#
        # Generating interpolation functions                              #
        #-----------------------------------------------------------------#
        ke_func     = Reginterp((x,y,z,time,), ke, method='linear') 
        enst_func   = Reginterp((x,y,z,time,), enst, method='linear') 
        #-----------------------------------------------------------------#
        # Extracting the values of interest                               #
        #-----------------------------------------------------------------#
        print('Extracting values of interest:')
        extract_values(ke_func, x1, x2, x3, t, data_path + 'ke.out')
        print('\tke-field')
        extract_values(enst_func, x1, x2, x3, t, data_path + 'enst.out')
        print('\tenst-field')
        #-----------------------------------------------------------------#
        # Deleting variables                                              #
        #-----------------------------------------------------------------#
        del ke
        del enst
        del ke_func
        del enst_func
    #---------------------------------------------------------------------#
    # CS2                                                                 #
    #---------------------------------------------------------------------#
    if cs2_flag is True:
        #-----------------------------------------------------------------#
        # Loading data                                                    #
        #-----------------------------------------------------------------#
        print('Load data:')
        time    = load(data_path + 'time.npy')[ts:tf]
        print('\ttime')
        cs2     = load(data_path + 'cs2.npy')[:,:,:,ts:tf]
        cs2     = ghost_cells(N, cs2)
        print('\tCS2')
        #-----------------------------------------------------------------#
        # Generating interpolation functions                              #
        #-----------------------------------------------------------------#
        cs2_func    = Reginterp((x,y,z,time,), cs2, method='linear') 
        #-----------------------------------------------------------------#
        # Extracting the values of interest                               #
        #-----------------------------------------------------------------#
        print('Extracting values of interest:')
        extract_values(cs2_func, x1, x2, x3, t, data_path + 'cs2.out')
        print('\tCS2')
        #-----------------------------------------------------------------#
        # Deleting variables                                              #
        #-----------------------------------------------------------------#
        del cs2
        del cs2_func

    print('**** Successful run ****')


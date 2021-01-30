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
from numpy import arange, mean, zeros, load, loadtxt, linspace, pi, flipud
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from ales_post.plot_settings    import plot_setting
from min_max                    import find_max3D
from rk4_back_tracking          import rk4_back_tracking
from rk4_back_tracking          import ghost_cells
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Generating dash                                                         #
#-------------------------------------------------------------------------#
def dash(
        number,
        length):
    """ Generating the line symbol """
    if number < len(length):
        line_not    = '*--'
    if number > len(length) and number < 2.0*len(length):
        line_not    = '--'
    else:
        line_not    = '^--'

    return line_not
#-------------------------------------------------------------------------#
# Generating label                                                        #
#-------------------------------------------------------------------------#
def get_label(
        val,
        switch):
    """ Generating the label for the location plots """
    #-----------------------------------------------------------------#
    # Generating label                                                #
    #-----------------------------------------------------------------#
    if int(val-switch) == 1:
        Lab = '$t_{0} + \\Delta t$' 
    elif int(val-switch) < 0:
        Lab = '$t_{0} - %i\\Delta t$'   %(abs(int(val-switch))) 
    elif int(val-switch) == 0:
        Lab = '$t_{0}$'    
    elif int(val-switch) > 1:
        Lab = '$t_{0} + %i\\Delta t$'   %(abs(int(val-switch))) 

    return Lab
#-------------------------------------------------------------------------#
# Tracking locations                                                      #
#-------------------------------------------------------------------------#
def plot_x_location(
        Locs,
        T_steps,
        T_vals,
        direction,
        name        = False,
        time_flag   = False):

    """ Plotting different x locations """

    #---------------------------------------------------------------------#
    # Default is time steps                                               #
    #---------------------------------------------------------------------#
    if time_flag is True:
        x_vec   = T_vals
        xlab    = 'Time'
    else:
        x_vec   = T_steps
        xlab    = 'Time step'
    #---------------------------------------------------------------------#
    # Plotting location                                                   #
    #---------------------------------------------------------------------#
    plot_setting()
    syms    = ['r', 'b', 'g', 'k', 'm', 'c', 'y']
    for q, T in enumerate(test):
        #-----------------------------------------------------------------#
        # Plot symbol                                                     #
        #-----------------------------------------------------------------#
        if q >= len(syms):
            diff    = int(q-len(syms))
            symbol  = syms[diff] + dash(diff, syms)
        else:
            symbol = syms[q]
        #-----------------------------------------------------------------#
        # Plot Label                                                      #
        #-----------------------------------------------------------------#
        lab     = get_label(T,tzero)
        #-----------------------------------------------------------------#
        # Plot plot                                                       #
        #-----------------------------------------------------------------#
        plt.plot(x_vec[q], Locs[q], symbol, label=lab)
        #-----------------------------------------------------------------#
        # Plotting the location of BS and DS switch                       #
        #-----------------------------------------------------------------#
        plt.plot([0, 0], [7,-7], 'b--', lw=1.5)
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
    #---------------------------------------------------------------------#
    plt.grid(True)
    plt.legend(loc=0)
    plt.ylabel('$x_{%i}$ location'          %(int(direction)))
    plt.xlabel(xlab)
    #---------------------------------------------------------------------#
    # Plot limits                                                         #
    #---------------------------------------------------------------------#
    x_min   = -int(T-tzero) 
    x_max   = int(max(test)-tzero) 
    plt.xlim([x_min, x_max])
    plt.ylim([0., 2.0*pi])
    plt.yticks([0, 0.5*pi, pi, 1.5*pi, 2.0*pi], \
                [0,\
                '$\\pi/2$',\
                '$\\pi$',\
                '$3\\pi/2$',\
                '$2\\pi$'])
    #---------------------------------------------------------------------#
    # Flip axis                                                           #
    #---------------------------------------------------------------------#
    ax1 = plt.gca()
    ax1.invert_xaxis()
    #---------------------------------------------------------------------#
    # Saving                                                              #
    #---------------------------------------------------------------------#
    if name is not False:
        plt.layout_tight()
        plt.savefig(media_path + name + '.png', bbox_inches='tight')
    plt.show()

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
    test        = [123, 124, 125]
    ts          = 0
    tf          = 150
    tzero       = 84
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
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    dx  = 2.*pi/N
    x1  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    x2  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    x3  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
    #---------------------------------------------------------------------#
    # Initializing strings                                                #
    #---------------------------------------------------------------------#
    x1_loc  = []
    x2_loc  = []
    x3_loc  = []
    t_val   = []
    t_steps = []
    #---------------------------------------------------------------------#
    # Test loop                                                           #
    #---------------------------------------------------------------------#
    IC      = '\n'
    for count,  tint in enumerate(test):
        #-----------------------------------------------------------------#
        # Maximum values and coordinates                                  #
        #-----------------------------------------------------------------#
        [val, xc, yc, zc]   = find_max3D(enst[:,:,:,tint])
        print('enst max --> %10.5e'                 %(val))
        print('x-index = %3i\tx-value = %8.4e'      %(xc, x[xc]))
        print('y-index = %3i\ty-value = %8.4e'      %(yc, y[yc]))
        print('z-index = %3i\tz-value = %8.4e'      %(zc, z[zc]))
        #-----------------------------------------------------------------#
        # Initial point                                                   #
        #-----------------------------------------------------------------#
        x0  = x[xc]
        y0  = y[yc]
        z0  = z[zc]
        t0  = time[tint]
        IC  += 'tval=%10i\t\tval=%20.8e\t\tx0=%20.8f\t\tx1=%20.8f\t\tx1=%20.8f\n'\
                %(tint, val, x0, y0, z0)
        #-----------------------------------------------------------------#
        # Domain variables                                                #
        #-----------------------------------------------------------------#
        dx  = 2.*pi/N
        x1  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
        x2  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
        x3  = linspace(-0.5*dx, 2.0*pi+0.5*dx, N+2)
        #-----------------------------------------------------------------#
        # Loading velocity                                                #
        #-----------------------------------------------------------------#
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
        #-----------------------------------------------------------------#
        # Initial start point                                             #
        #-----------------------------------------------------------------#
        print('%25.20f'     %(t0))
        print('%25.20f'     %(x0))
        print('%25.20f'     %(y0))
        print('%25.20f'     %(z0))
        #-----------------------------------------------------------------#
        # Back tracking                                                   #
        #-----------------------------------------------------------------#
        rk4_back_tracking([x0, y0, z0, t0], u1, u2, u3, \
                            x1, x2, x3, time, tint, 0,\
                            'enstrophy-back-print.out', \
                            'enstrophy-back-coordinates.out')
        #-----------------------------------------------------------------#
        # Storing x1-x3                                                   #
        #-----------------------------------------------------------------#
        data    = loadtxt('enstrophy-back-coordinates.out',skiprows=1)
        x1_loc.append(data[:,1])
        x2_loc.append(data[:,2])
        x3_loc.append(data[:,3])
        t_val.append(data[:,0])
        t_steps.append(flipud(linspace(ts-tzero, tint-tzero, len(x1_loc[count]))))
    #---------------------------------------------------------------------#
    # Printing initial conditions                                         #
    #---------------------------------------------------------------------#
    call(['clear'])
    print(IC)
    #---------------------------------------------------------------------#
    # Plotting data                                                       #
    #---------------------------------------------------------------------#
    plot_x_location( x1_loc, t_steps, t_val, 1)
    plot_x_location( x2_loc, t_steps, t_val, 2)
    plot_x_location( x3_loc, t_steps, t_val, 3)
    #---------------------------------------------------------------------#
    # Printing initial conditions                                         #
    #---------------------------------------------------------------------#
    call(['clear'])
    print(IC)

    print('**** Successful run ****')
    sys.exit(0)

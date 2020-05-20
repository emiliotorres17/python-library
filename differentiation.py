#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate gradients using either
    spectral or finite difference tools.

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
import matplotlib.pyplot as plt 
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Spectral differentiation                                                #
#-------------------------------------------------------------------------#
def spec_diff(
        field,
        grad_dir):

    """ Calculating gradients using spectral methods """
    #---------------------------------------------------------------------#
    # Setting the k-space field                                           #
    #---------------------------------------------------------------------#
    dim     = field.shape[0] 
    k       = np.fft.fftfreq(dim) * dim
    Kfield  = np.array(np.meshgrid(k,k,k,indexing='ij'))
    #---------------------------------------------------------------------#
    # Differentiating                                                     #
    #---------------------------------------------------------------------#
    derv    = np.fft.ifftn(1j*Kfield[grad_dir]*np.fft.fftn(field)).real

    return derv
#-------------------------------------------------------------------------#
# Gradient differentiation                                                #
#-------------------------------------------------------------------------#
def grad_diff(      
        field,                          # field to differentiate
        grad_dir,                       # gradient direction
        step_size):                     # spatial step size

    """ Calculating gradients using finite differencing """ 
    #---------------------------------------------------------------------#
    # Differentiating                                                     #
    #---------------------------------------------------------------------#
    derv    = np.gradient( field, step_size, edge_order=2)[grad_dir] 

    return derv
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    print('unit test is visual')
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 200
    x           = np.linspace(0.0, 2.0*np.pi, N)
    y           = np.linspace(0.0, 2.0*np.pi, N)
    z           = np.linspace(0.0, 2.0*np.pi, N)
    dx          = (2.0*np.pi)/(N-1)
    [X1, X2]    = np.meshgrid(x,x,)
    f           = np.zeros((N,N,N))
    dfdx        = np.zeros((N,N,N))
    dfdy        = np.zeros((N,N,N))
    dfdz        = np.zeros((N,N,N))
    #---------------------------------------------------------------------#
    # Plotting solutions                                                  #
    #---------------------------------------------------------------------#
    c1          = 1.0
    vmin        = -c1
    vmax        = c1
    dp          = (vmax-vmin)/500.00
    ticks       = [0,3,6]
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    count = 0
    for k in range(0, N):
        for j in range(0,N):
            for i in range(0,N):
                f[i,j,k]    = np.sin(x[i])*np.sin(y[j])*np.sin(z[k])
                dfdx[i,j,k] = np.cos(x[i])*np.sin(y[j])*np.sin(z[k])    
                dfdy[i,j,k] = np.sin(x[i])*np.cos(y[j])*np.sin(z[k])    
                dfdz[i,j,k] = np.sin(x[i])*np.sin(y[j])*np.cos(z[k])    
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if count > 20:
            print('exact solution --> %i'       %(k))
            count   = 0
        count += 1
    #---------------------------------------------------------------------#
    # Spectral solution                                                   #
    #---------------------------------------------------------------------#
    spec_dx     = spec_diff(f, 0)
    spec_dy     = spec_diff(f, 1)
    spec_dz     = spec_diff(f, 2)
    print('Calculated spectral solution')
    #---------------------------------------------------------------------#
    # Gradient solution                                                   #
    #---------------------------------------------------------------------#
    grad_dx     = grad_diff(f, 0, dx)
    grad_dy     = grad_diff(f, 1, dx)
    grad_dz     = grad_diff(f, 2, dx)
    print('Calculated gradient solution')
    #---------------------------------------------------------------------#
    # Plotting solution                                                   #
    #---------------------------------------------------------------------#
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    ax  = plt.subplot(1,3,1, aspect='equal')
    cnt = plt.contourf(X1, X2, dfdx[:,:,31], np.arange(vmin, vmax, dp),\
                    cmap='jet', extend='both')
    ax.set_yticks(ticks)
    ax.set_xticks(ticks)
    for c in cnt.collections:
        c.set_edgecolor('face')
    print('Plotted exact solution')
    #---------------------------------------------------------------------#
    #  Spectral solution                                                  #
    #---------------------------------------------------------------------#
    ax  = plt.subplot(1,3,2, aspect='equal')
    cnt = plt.contourf(X1, X2, spec_dx[:,:,31], np.arange(vmin, vmax, dp),\
                    cmap='jet', extend='both')
    ax.set_yticks(ticks)
    ax.set_xticks(ticks)
    for c in cnt.collections:
        c.set_edgecolor('face')
    print('Plotted spectral solution')
    #---------------------------------------------------------------------#
    #  Gradient solution                                                  #
    #---------------------------------------------------------------------#
    ax  = plt.subplot(1,3,3, aspect='equal')
    cnt = plt.contourf(X1, X2, grad_dx[:,:,31], np.arange(vmin, vmax, dp),\
                    cmap='jet', extend='both')
    ax.set_yticks(ticks)
    ax.set_xticks(ticks)
    for c in cnt.collections:
        c.set_edgecolor('face')
    print('Plotted gradient solution')
    
    cax = plt.axes([0.925, 0.3467, 0.0125, 0.30])
    plt.colorbar(cax=cax, ticks=[-1, 0, 1])
    plt.show()
    
    print('**** Successful run ****')
    sys.exit(0)

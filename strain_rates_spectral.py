#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the strain rates using
    differentiation.

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
# Spectral strain rates                                                   #
#-------------------------------------------------------------------------#
def spectral_strain_rates(
        U1,                     # velocity-1
        U2,                     # velocity-2
        U3):                    # velocity-3

    """ Calculating the strain rates using spectral methods """
    #---------------------------------------------------------------------#
    # Defining the domain variables                                       #
    #---------------------------------------------------------------------#
    dim     = U1.shape[0]
    kspec   = np.fft.fftfreq(dim) * dim
    Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    St      = np.zeros((6, dim, dim, dim))
    U       = np.zeros((3, dim, dim, dim))
    U[0]    = U1
    U[1]    = U2
    U[2]    = U3
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    St[0]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[0]) + 1j*Kfield[0]*np.fft.fftn(U[0])).real
    St[1]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[1]) + 1j*Kfield[1]*np.fft.fftn(U[0])).real
    St[2]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[0])).real
    St[3]   = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U[1]) + 1j*Kfield[1]*np.fft.fftn(U[1])).real
    St[4]   = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[1])).real
    St[5]   = 0.5*np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(U[2]) + 1j*Kfield[2]*np.fft.fftn(U[2])).real
    #St[0]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U1) + 1j*Kfield[0]*np.fft.fftn(U1)).real
    #St[1]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U2) + 1j*Kfield[1]*np.fft.fftn(U1)).real
    #St[2]   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(U3) + 1j*Kfield[2]*np.fft.fftn(U1)).real
    #St[3]   = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U2) + 1j*Kfield[1]*np.fft.fftn(U2)).real
    #St[4]   = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(U3) + 1j*Kfield[2]*np.fft.fftn(U2)).real
    #St[5]   = 0.5*np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(U3) + 1j*Kfield[2]*np.fft.fftn(U3)).real

    return St
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    pi      = np.pi
    N       = 32
    x0      = 0.0
    xf      = 2.0*pi
    x       = np.linspace(x0, xf, N+1)
    y       = np.linspace(x0, xf, N+1)
    z       = np.linspace(x0, xf, N+1)
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    u       = np.zeros((N+1, N+1, N+1))
    v       = np.zeros((N+1, N+1, N+1))
    w       = np.zeros((N+1, N+1, N+1))
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Calculating velocities                                              #
    #---------------------------------------------------------------------#
    print_count = 0
    print("**** Velocity loop")
    for k in range(0,N+1):
        for j in range(0,N+1):
            for i in range(0,N+1):
                u[i,j,k]            = np.cos(pi*x[i])
                v[i,j,k]            = x[i]*z[k]*np.sin(pi*y[j])
                w[i,j,k]            = x[i]*y[j]*np.cos(pi*z[k])
        #-----------------------------------------------------------------#
        # Print statement output                                          #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('z = %i'          %(k))
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate solution                                #
    #---------------------------------------------------------------------#
    S   = spectral_strain_rates(u, v, w)[0]
    print("**** Calculated approximate strain rate ****")
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del u
    del v
    del w
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating the strain rates                                      #
    #---------------------------------------------------------------------#
    S11_exact   = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating strain rates                                            #
    #---------------------------------------------------------------------#
    print("**** Strain rate loop:")
    for k in range(0,N+1):
        for j in range(0,N+1):
            for i in range(0,N+1):
                S11_exact[i,j,k]    = -pi*np.sin(pi*x[i])
        #-----------------------------------------------------------------#
        # Print statement output                                          #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('z = %i'              %(k))
            print_count = 0
        print_count += 1
    #=====================================================================#
    # Post processing                                                     #
    #=====================================================================#
    zslice      = int((N+1)/2)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Dummy plot                                                          #
    #---------------------------------------------------------------------#
    dummy       = np.zeros((N+1, N+1))
    dummy[0,0]  = np.amax(S11_exact[:,:,zslice])
    print(np.amax(S11_exact[:,:,zslice]))
    dummy[0,1]  = np.amin(S11_exact[:,:,zslice])
    print(np.amin(S11_exact[:,:,zslice]))
    cnt1        = plt.contourf(X1, X2, dummy, 500, cmap="jet") 
    plt.colorbar()
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    cnt         = plt.contourf(X1, X2, S[:,:,zslice], 500, cmap="jet") 
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0\leq x \leq \pi$")
    plt.ylabel("$0\leq y \leq \pi$")
    plt.colorbar(cnt1)
    plt.savefig("approximate.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt         = plt.contourf(X1, X2, S11_exact[:,:,zslice], 500, cmap="jet") 
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0\leq x \leq \pi$")
    plt.ylabel("$0\leq y \leq \pi$")
    plt.colorbar(cnt1)
    plt.savefig("exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error plot                                                          #
    #---------------------------------------------------------------------#
    error       = abs(S - S11_exact)
    cnt         = plt.contourf(X1[20:30,20:30], X2[20:30,20:30],\
                        error[20:30,20:30,zslice], 500, cmap="jet") 
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0\leq x \leq \pi$")
    plt.ylabel("$0\leq y \leq \pi$")
    plt.colorbar()
    plt.savefig("error.pdf")
    plt.clf()

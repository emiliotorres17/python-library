#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate the A term, vortex
    stretching term, in the enstrophy transport term using both the
    spectral and numpy gradient tools.

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
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from strain_rates_spectral      import spectral_strain_rates
from strain_rates               import strain_rates
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# A term (Vortex stretching)                                              #
#-------------------------------------------------------------------------#
def a_term_enstrophy(
        U1,                     # velocity 1 component
        U2,                     # velocity 2 component
        U3,                     # velocity 3 component
        W1,                     # vorticity 1 component
        W2,                     # vorticity 2 component
        W3,                     # vorticity 3 component
        flag    = False,        # spectral flag
        h       = False,        # spatial step size
        dim     = False):       # Number of spatial steps
    """ Calculating the vortex stretching term in enstrophy transport
        equation """
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    if h is False:
        Pi  = np.pi
        num = 64
        h   = (2.0*Pi)/num
    if dim is False:
        dim = 64
    a   = np.zeros((dim,dim,dim))
    #---------------------------------------------------------------------#
    # Strain rates                                                        #
    #---------------------------------------------------------------------#
    if flag is False:
        St  = spectral_strain_rates(U1, U2, U3)
    else:
        St  = strain_rates(U1,U2,U3, h, dim)
    #---------------------------------------------------------------------#
    # Calculating vortex stretching                                       #
    #---------------------------------------------------------------------#
    a   += W1*St[0]*W1
    a   += 2.0*W1*St[1]*W2
    a   += 2.0*W3*St[2]*W1
    a   += W2*St[3]*W2
    a   += 2.0*W2*St[4]*W3
    a   += W3*St[5]*W3

    return a
#=========================================================================#
# Main preamble                                                           #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main                                                                #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep             = os.sep
    pwd             = os.getcwd()
    media_path      = pwd + "%cmedia%c"             %(sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    N           = 400
    dx          = (2.0*pi)/(N)
    x           = np.linspace(0, 2.0*pi, N)
    y           = np.linspace(0, 2.0*pi, N)
    z           = np.linspace(0, 2.0*pi, N)
    (X1, X2)    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    u1  = np.zeros((N,N,N))
    u2  = np.zeros((N,N,N))
    u3  = np.zeros((N,N,N))
    w1  = np.zeros((N,N,N))
    w2  = np.zeros((N,N,N))
    w3  = np.zeros((N,N,N))
    for k in range(0,N):
        for j in range(0,N):
            for i in range(0,N):
                u1[i,j,k]   = np.cos(x[i])
                u2[i,j,k]   = np.sin(y[j])
                u3[i,j,k]   = np.cos(z[k])
                w1[i,j,k]   = np.sin(x[i])
                w2[i,j,k]   = np.cos(y[j])
                w3[i,j,k]   = np.sin(z[k])
        print(k)
    A       = a_term_enstrophy(u1, u2, u3, w1, w2, w3, False, dx, N)
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del u1
    del u2
    del u3
    del w1
    del w2
    del w3
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    sol = np.zeros((N,N,N))
    for k in range(0,N):
        for j in range(0,N):
            for i in range(0,N):
                sol[i,j,k] = -(np.sin(x[i]))**3.0 + (np.cos(y[j]))**3.0 -\
                                (np.sin(z[k]))**3.0
        print(k)
    #---------------------------------------------------------------------#
    # Plotting exact solution                                             #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1,X2, sol[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting approximate solution                                       #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1,X2, A[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting the error                                                  #
    #---------------------------------------------------------------------#
    error   = abs(sol-A)
    cnt     = plt.contourf(X1,X2, error[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "A-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

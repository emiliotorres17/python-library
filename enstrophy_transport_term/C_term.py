#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script to calculate the C term, laplacian(Omega),
    in the enstrophy transport equation.

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
# C term                                                                  #
#-------------------------------------------------------------------------#
def c_term_enstrophy(
        enst,                           # enstrophy field
        h       = False,                # spatial step size
        Nu      = False,                # viscosity
        flag    = True):                # spectral flag; default is gradient

    """ Calculating the diffusion term in enstrophy transport equation """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if h is False:
        num = 64
        Pi  = np.pi
        h  = (2.0*pi)/num
    if Nu is False:
        Nu  = 0.000185
    #---------------------------------------------------------------------#
    # Calculating the Laplacian                                           #
    #---------------------------------------------------------------------#
    dim = enst.shape[0]
    if flag is True:
        term1   = np.gradient(np.gradient(enst, h, edge_order=2)[0],\
                                h, edge_order=2)[0]
        term2   = np.gradient(np.gradient(enst, h, edge_order=2)[1],\
                                h, edge_order=2)[1]
        term3   = np.gradient(np.gradient(enst, h, edge_order=2)[2],\
                                h, edge_order=2)[2]
    else:
        term1   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term1   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term1)).real
        term2   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term2   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term2)).real
        term3   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst)).real
        term3   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(term3)).real
    #---------------------------------------------------------------------#
    # Calculating dissipation term                                        #
    #---------------------------------------------------------------------#
    c   = np.zeros((dim,dim,dim))
    c   += Nu*(term1 + term2 + term3)

    return c
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep         = os.sep
    pwd         = os.getcwd()
    media_path  = pwd + "%cmedia%c"         %(sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    N           = 128
    x0          = -1.0
    xf          = 1.0
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    enst        = np.zeros((N+1,N+1,N+1))
    sol         = np.zeros((N+1,N+1,N+1))
    #---------------------------------------------------------------------#
    # Calculating approximate and exact solution                          #
    #---------------------------------------------------------------------#
    for k in range(0,N+1):
        for j in range(0,N+1):
            for i in range(0,N+1):
                enst[i,j,k] = x[i]**2.0*y[j]**2.0*z[k]**2.0
                sol[i,j,k]  = 2.0*y[j]**2.0*z[k]**2.0 +\
                                2.0*x[i]**2.0*z[k]**2.0 +\
                                2.0*x[i]**2.0*y[j]**2.0
        print(k)
    c_approx    = c_term_enstrophy(enst, dx, 1.0, True)
    #---------------------------------------------------------------------#
    # Plotting approximate solution                                       #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, c_approx[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting exact solution                                             #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, sol[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error                                                               #
    #---------------------------------------------------------------------#
    error   = abs(sol-c_approx)
    cnt     = plt.contourf(X1, X2, error[:,:,31], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

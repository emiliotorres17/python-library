#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the dissipation term
    (D term) in he enstrophy transport term.

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
def d_term_enstrophy(
        w1,                 # vorticity component-1
        w2,                 # vorticity component-2
        w3,                 # vorticity component-3
        h       = False,    # spatial step size
        Nu      = False,    # viscosity
        flag    = True):    # spectral flag; default is gradient tool

    """ Calculating the enstrophy dissipation term """
    #---------------------------------------------------------------------#
    # Default settings                                                    #
    #---------------------------------------------------------------------#
    if Nu is False:             # default viscosity
        Nu  = 0.000185
    if h is False:
        Pi  = np.pi
        num = 64
        h   = (2.0*Pi)/num
    #---------------------------------------------------------------------#
    # Calculating the gradients                                           #
    #---------------------------------------------------------------------#
    dim     = w1.shape[0]
    if flag is True:
        term1   = np.gradient(w1, h, edge_order=2)[0]
        term2   = np.gradient(w1, h, edge_order=2)[1]
        term3   = np.gradient(w1, h, edge_order=2)[2]
        term4   = np.gradient(w2, h, edge_order=2)[0]
        term5   = np.gradient(w2, h, edge_order=2)[1]
        term6   = np.gradient(w2, h, edge_order=2)[2]
        term7   = np.gradient(w3, h, edge_order=2)[0]
        term8   = np.gradient(w3, h, edge_order=2)[1]
        term9   = np.gradient(w3, h, edge_order=2)[2]
    else:
        kspec   = np.fft.fftfreq(dim) * dim
        Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
        term1   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w1)).real
        term2   = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w1)).real
        term3   = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w1)).real
        term4   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w2)).real
        term5   = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w2)).real
        term6   = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w2)).real
        term7   = np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w3)).real
        term8   = np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w3)).real
        term9   = np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w3)).real
    #---------------------------------------------------------------------#
    # Calculating the dissipation                                         #
    #---------------------------------------------------------------------#
    d   = np.zeros((dim, dim, dim))
    d   += -Nu*(term1*term1 + term2*term2 + term3*term3 +\
                term4*term4 + term5*term5 + term6*term6 +\
                term7*term7 + term8*term8 + term9*term9)

    return d
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    sep             = os.sep
    pwd             = os.getcwd()
    media_path      = pwd + "%cmedia%c"         %(sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    nu          = 1.0
    N           = 400
    pi          = np.pi
    dx          = (2.0*pi)/N
    x           = np.linspace(0, 2.0*pi, N)
    y           = np.linspace(0, 2.0*pi, N)
    z           = np.linspace(0, 2.0*pi, N)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Pre-allocating space                                                #
    #---------------------------------------------------------------------#
    omega1  = np.zeros((N,N,N))
    omega2  = np.zeros((N,N,N))
    omega3  = np.zeros((N,N,N))
    sol     = np.zeros((N,N,N))
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    for k in range(0,N):
        for j in range(0,N):
            for i in range(0,N):
                omega1[i,j,k]   = np.sin(x[i])
                omega2[i,j,k]   = np.cos(y[j])
                omega3[i,j,k]   = np.sin(z[k])
                sol[i,j,k]      = -nu*(np.cos(x[i])**2.0 +\
                                    np.sin(y[j])**2.0 +\
                                    np.cos(z[k])**2.0)
        print(k)
    d_approx = d_term_enstrophy(omega1, omega2, omega3, dx, nu, True)
    #---------------------------------------------------------------------#
    # Plotting approximate solution                                       #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, d_approx[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "D-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting exact solution                                             #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, sol[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "D-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting error                                                      #
    #---------------------------------------------------------------------#
    error   = abs(sol-d_approx)
    cnt     = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "D-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

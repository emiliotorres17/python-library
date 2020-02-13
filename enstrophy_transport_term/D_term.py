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
    d       = np.zeros((dim, dim, dim))
    if flag is True:
        grad1   = np.gradient(w1, h, edge_order=2)
        grad2   = np.gradient(w2, h, edge_order=2)
        grad3   = np.gradient(w3, h, edge_order=2)
        
        d       += (grad1[0])**2.0
        d       += (grad1[1])**2.0
        d       += (grad1[2])**2.0
        d       += (grad2[0])**2.0
        d       += (grad2[1])**2.0
        d       += (grad2[2])**2.0
        d       += (grad3[0])**2.0
        d       += (grad3[1])**2.0
        d       += (grad3[2])**2.0
    else:
        kspec   = np.fft.fftfreq(dim) * dim
        Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
        d       += (np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w1)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w1)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w1)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w2)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w2)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w2)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(w3)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(w3)).real)**2.0
        d       += (np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(w3)).real)**2.0
    #---------------------------------------------------------------------#
    # Calculating the dissipation                                         #
    #---------------------------------------------------------------------#
    d   *= -Nu

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
    N           = 256
    xf          = 1.0
    x0          = 0.0
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Pre-allocating space                                                #
    #---------------------------------------------------------------------#
    omega1  = np.zeros((N+1,N+1,N+1))
    omega2  = np.zeros((N+1,N+1,N+1))
    omega3  = np.zeros((N+1,N+1,N+1))
    sol     = np.zeros((N+1,N+1,N+1))
    #---------------------------------------------------------------------#
    # Approximate and exact solutions                                     #
    #---------------------------------------------------------------------#
    for k in range(0,N+1):
        for j in range(0,N+1):
            for i in range(0,N+1):
                #---------------------------------------------------------#
                # Vorticities                                             #
                #---------------------------------------------------------#
                omega1[i,j,k]   = y[j]*z[k]*np.sin(x[i])
                omega2[i,j,k]   = x[i]*z[k]*np.sin(y[j])
                omega3[i,j,k]   = x[i]*y[j]*np.sin(z[k])
                sol[i,j,k]      = (y[j]*z[k]*np.cos(x[i]))**2.0\
                                    + (z[k]*np.sin(x[i]))**2.0\
                                    + (y[j]*np.sin(x[i]))**2.0\
                                    + (z[k]*np.sin(y[j]))**2.0\
                                    + (x[i]*z[k]*np.cos(y[j]))**2.0\
                                    + (x[i]*np.sin(y[j]))**2.0\
                                    + (y[j]*np.sin(z[k]))**2.0\
                                    + (x[i]*np.sin(z[k]))**2.0\
                                    + (y[j]*x[i]*np.cos(z[k]))**2.0
        print(k)
    d_approx = d_term_enstrophy(omega1, omega2, omega3, dx, nu, True)
    sol *= -nu
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

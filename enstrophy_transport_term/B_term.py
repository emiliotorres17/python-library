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
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from enstrophy      import enstrophy_static
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# B term                                                                  #
#-------------------------------------------------------------------------#
def b_term_enstrophy(
        w1,                             # vorticity-1 component
        w2,                             # vorticity-2 component
        w3,                             # vorticity-3 component
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
        h  = (2.0*Pi)/num
    if Nu is False:
        Nu  = 0.000185
    #---------------------------------------------------------------------#
    # Calculating the enstrophy                                           #
    #---------------------------------------------------------------------#
    Enst    = enstrophy_static(w1, w2, w3)
    #---------------------------------------------------------------------#
    # Calculating the Laplacian                                           #
    #---------------------------------------------------------------------#
    dim = Enst.shape[0]
    b   = np.zeros((dim,dim,dim))
    if flag is True:
        b       += np.gradient(np.gradient(Enst, h, edge_order=2)[0],\
                                h, edge_order=2)[0]
        b       += np.gradient(np.gradient(Enst, h, edge_order=2)[1],\
                                h, edge_order=2)[1]
        b       += np.gradient(np.gradient(Enst, h, edge_order=2)[2],\
                                h, edge_order=2)[2]
        b       *= Nu
    else:
        kspec   = np.fft.fftfreq(dim) * dim
        Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
        Ksq     = np.sum(np.square(Kfield), axis=0)
        term1   = -Nu*Ksq*np.fft.fftn(Enst)
        b       = np.fft.ifftn(term1).real

    return b
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
    N           = 256
    x0          = 0
    xf          = 1.0
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    omega1      = np.zeros((N+1,N+1,N+1))
    omega2      = np.zeros((N+1,N+1,N+1))
    omega3      = np.zeros((N+1,N+1,N+1))
    sol         = np.zeros((N+1,N+1,N+1))
    #---------------------------------------------------------------------#
    # Calculating approximate and exact solution                          #
    #---------------------------------------------------------------------#
    for k in range(0,N+1):
        for j in range(0,N+1):
            for i in range(0,N+1):
                #---------------------------------------------------------#
                # vorticities                                             #
                #---------------------------------------------------------#
                omega1[i,j,k]   = y[j]*z[k]*np.sin(x[i])
                omega2[i,j,k]   = x[i]*z[k]*np.sin(y[j])
                omega3[i,j,k]   = x[i]*y[j]*np.sin(z[k])
                #---------------------------------------------------------#
                # exact solution                                          #
                #---------------------------------------------------------#
                df_dx       = y[j]**2.0*z[k]**2.0*\
                                (4.0*np.cos(x[i])**2.0-2.0)\
                                + 2.0*z[k]**2.0*np.sin(y[j])**2.0\
                                + 2.0*y[j]**2.0*np.sin(z[k])**2.0
                df_dy       = x[i]**2.0*z[k]**2.0*\
                                (4.0*np.cos(y[j])**2.0-2.0)\
                                + 2.0*z[k]**2.0*np.sin(x[i])**2.0\
                                + 2.0*x[i]**2.0*np.sin(z[k])**2.0
                df_dz       = x[i]**2.0*y[j]**2.0*\
                                (4.0*np.cos(z[k])**2.0-2.0)\
                                + 2.0*y[j]**2.0*np.sin(x[i])**2.0\
                                + 2.0*x[i]**2.0*np.sin(y[j])**2.0
                sol[i,j,k]  = 0.5*(df_dx + df_dy + df_dz)
        print(k)
    c_approx    = c_term_enstrophy(omega1, omega2, omega3, dx, 1.0, True)
    #---------------------------------------------------------------------#
    # Plotting approximate solution                                       #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, c_approx[:,:,31], 500, cmap='jet')
    for C in cnt.collections:
        C.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting exact solution                                             #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, sol[:,:,31], 500, cmap='jet')
    for C in cnt.collections:
        C.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error                                                               #
    #---------------------------------------------------------------------#
    error   = abs(sol-c_approx)
    cnt     = plt.contourf(X1, X2, error[:,:,31], 500, cmap='jet')
    for C in cnt.collections:
        C.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "C-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

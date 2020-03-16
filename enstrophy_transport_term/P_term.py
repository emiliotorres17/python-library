#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the SGS production in
    the enstrophy transport equation.

    **** Note: The spectral differentiation still needs to be added

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
from enstrophy_transport_term.psi_enstrophy     import psi
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
def p_term_enstrophy(
        w1,                 # vorticity-1 component
        w2,                 # vorticity-2 component
        w3,                 # vorticity-3 component
        Tau,                # SGS; (6,64,64,64)
        h       = False,    # spatial step size
        flag    = True):    # spectral flag; default gradient tool

    """ Calculating the SGS enstrophy production """
    #---------------------------------------------------------------------#
    # Default variables                                                   #
    #---------------------------------------------------------------------#
    if h is False:
        Pi  = np.pi
        num = 64
        h   = (2.0*Pi)/num
    #---------------------------------------------------------------------#
    # Calculating psi                                                     #
    #---------------------------------------------------------------------#
    Psi = psi(Tau, h, flag)
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    dim         = np.shape(w1)[0]
    SGS_prod    = np.zeros((dim, dim, dim))
    #---------------------------------------------------------------------#
    # Calculating gradients                                               #
    #---------------------------------------------------------------------#
    grad_w1     = np.gradient(w1, h, edge_order=2)
    grad_w2     = np.gradient(w2, h, edge_order=2)
    grad_w3     = np.gradient(w3, h, edge_order=2)
    #---------------------------------------------------------------------#
    # Calculating SGS production                                          #
    #---------------------------------------------------------------------#
    if flag is False:               # spectral method
        print("Spectral differentiation has not been added")
        sys.exit(1)
    else:                           # gradient tool method
        #-----------------------------------------------------------------#
        # terms 1-3 (i=1)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[0]*grad_w1[0]
        SGS_prod    += Psi[1]*grad_w1[1]
        SGS_prod    += Psi[2]*grad_w1[2]
        #-----------------------------------------------------------------#
        # terms 4-6 (i=2)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[3]*grad_w2[0]
        SGS_prod    += Psi[4]*grad_w2[1]
        SGS_prod    += Psi[5]*grad_w2[2]
        #-----------------------------------------------------------------#
        # terms 7-9 (i=3)                                                 #
        #-----------------------------------------------------------------#
        SGS_prod    += Psi[6]*grad_w3[0]
        SGS_prod    += Psi[7]*grad_w3[1]
        SGS_prod    += Psi[8]*grad_w3[2]

    return SGS_prod
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
    N           = 128
    pi          = np.pi
    x0          = 0.0
    xf          = 1.0
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x,y)
    #---------------------------------------------------------------------#
    # Preallocating space                                                 #
    #---------------------------------------------------------------------#
    sol         = np.zeros((N+1, N+1, N+1))
    omega1      = np.zeros((N+1, N+1, N+1))
    omega2      = np.zeros((N+1, N+1, N+1))
    omega3      = np.zeros((N+1, N+1, N+1))
    tau         = np.zeros((6,N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Looping over domain                                                 #
    #---------------------------------------------------------------------#
    count   = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                #---------------------------------------------------------#
                # Vorticities                                             #
                #---------------------------------------------------------#
                omega1[i,j,k]   = y[j]*z[k]*np.sin(x[i])
                omega2[i,j,k]   = x[i]*z[k]*np.sin(y[j])
                omega3[i,j,k]   = x[i]*y[j]*np.sin(z[k])
                #---------------------------------------------------------#
                # SGS                                                     #
                #---------------------------------------------------------#
                tau[0,i,j,k]    = x[i]**2.0
                tau[1,i,j,k]    = x[i]*y[j]
                tau[2,i,j,k]    = x[i]*z[k]
                tau[3,i,j,k]    = y[j]**2.0
                tau[4,i,j,k]    = y[j]*z[k]
                tau[5,i,j,k]    = z[k]**2.0
                #---------------------------------------------------------#
                # Exact solution                                          #
                #---------------------------------------------------------#
                sol[i,j,k]      =   z[k]**2.0*np.sin(x[i])\
                                  - y[j]**2.0*np.sin(x[i])\
                                  - z[k]**2.0*np.sin(y[j])\
                                  + x[i]**2.0*np.sin(y[j])\
                                  + y[j]**2.0*np.sin(z[k])\
                                  - x[i]**2.0*np.sin(z[k])
        #-----------------------------------------------------------------#
        # Print criteria                                                  #
        #-----------------------------------------------------------------#
        if count == 20:
            print(k)
            count = 0
        count += 1
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    sgs_approx  = p_term_enstrophy(omega1, omega2, omega3, tau,\
                    dx, True)
    #---------------------------------------------------------------------#
    # Plotting approximate solution                                       #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, sgs_approx[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "SGS-prod-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting exact solution                                             #
    #---------------------------------------------------------------------#
    cnt     = plt.contourf(X1, X2, sol[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "SGS-prod-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Plotting error                                                      #
    #---------------------------------------------------------------------#
    error   = abs(sol-sgs_approx)
    cnt     = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar()
    plt.savefig(media_path + "SGS-prod-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the right hand side of
    the kinetic energy transport equation.

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
from central_difference         import time_derv
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Subroutine to calculate the RHS for the KE transport                    #
#-------------------------------------------------------------------------#
def KE_RHS(
        ke,                     # kinetic energy field --> NxNxNxN_time
        u1,                     # velocity-1 field --> NxNxN
        u2,                     # velocity-2 field --> NxNxN
        u3,                     # velocity-3 field --> NxNxN
        time,                   # time vector --> N_time
        Nt,                     # time step of interest
        h):                     # spatial step size

    """ Calculating the RHS of the kinetic energy transport equation """
    #---------------------------------------------------------------------#
    # Calculating the time derivative                                     #
    #---------------------------------------------------------------------#
    dt_term     = time_derv(ke, Nt, time)
    #delta_t     = time[Nt] - time[Nt-1]
    #dt_term     = np.gradient(ke, delta_t, edge_order=2)[3]
    #dt_term     = dt_term[:,:,:,Nt]
    #print(dt_term.shape)
    #sys.exit()
    #---------------------------------------------------------------------#
    # Calculating the advection terms                                     #
    #---------------------------------------------------------------------#
    term1       = u1*np.gradient(ke[:,:,:,Nt], h, edge_order=2)[0]
    term2       = u2*np.gradient(ke[:,:,:,Nt], h, edge_order=2)[1]
    term3       = u3*np.gradient(ke[:,:,:,Nt], h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating the RHS                                                 #
    #---------------------------------------------------------------------#
    RHS         = dt_term # + term1 + term2 + term3

    return RHS
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(["clear"])
    sep                 = os.sep
    pwd                 = os.getcwd()
    media_path          = pwd + "%cmedia%c"         %(sep, sep)
    #=====================================================================#
    # Defining domain variables                                           #
    #=====================================================================#
    pi      = np.pi
    N       = 64
    num_t   = 1500
    x0      = 0.0
    xf      = 1.0
    dx      = (xf-x0)/N
    t0      = 0.0
    tf      = 200.0
    t_int   = int(num_t - 4)
    x       = np.linspace(x0, xf, N+1)
    y       = np.linspace(x0, xf, N+1)
    z       = np.linspace(x0, xf, N+1)
    t       = np.linspace(t0, tf, num_t)
    #---------------------------------------------------------------------#
    # Preallocating the kinetic energy, ux, uy, uz                        #
    #---------------------------------------------------------------------#
    ke      = np.zeros((N+1, N+1, N+1, num_t))
    ux      = np.zeros((N+1, N+1, N+1, num_t))
    uy      = np.zeros((N+1, N+1, N+1, num_t))
    uz      = np.zeros((N+1, N+1, N+1, num_t))
    #=====================================================================#
    # Calculating the approximate Dk/Dt                                   #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Calculating the kinetic energy, ux, uy, uz                          #
    #---------------------------------------------------------------------#
    print_count = 0
    for n in range(0, num_t):
        for k in range(0, N+1):
            for j in range(0, N+1):
                for i in range(0, N+1):
                    ke[i,j,k,n]     = x[i]**2.0 * y[j]**2.0 * z[k]**2.0\
                                        * np.cos(2.0*pi*t[n])
                    ux[i,j,k,n]     = x[i]*np.sin(pi*t[n])
                    uy[i,j,k,n]     = y[j]*np.cos(pi*t[n])
                    uz[i,j,k,n]     = z[k]*np.sin(pi*t[n])
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print("Approximation loop --> t_step = %i"      %(n))
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate RHS                                     #
    #---------------------------------------------------------------------#
    RHS_approx      = KE_RHS(ke, ux[:,:,:,t_int], uy[:,:,:,t_int],\
                                uz[:,:,:,t_int], t, t_int, dx)
    print("**** Calculated the approximate RHS ****")
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del ke
    del ux
    del uy
    del uz
    #=====================================================================#
    # Calculating the exact solution                                      #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating terms for the exact solution                          #
    #---------------------------------------------------------------------#
    term_t  = np.zeros((N+1, N+1, N+1, num_t))
    term_x  = np.zeros((N+1, N+1, N+1, num_t))
    term_y  = np.zeros((N+1, N+1, N+1, num_t))
    term_z  = np.zeros((N+1, N+1, N+1, num_t))
    #---------------------------------------------------------------------#
    # Calculating the kinetic energy, ux, uy, uz                          #
    #---------------------------------------------------------------------#
    print_count = 0
    for n in range(0, num_t):
        for k in range(0, N+1):
            for j in range(0, N+1):
                for i in range(0, N+1):
                    term_t[i,j,k,n]     = -2.0*pi*x[i]**2.0*y[j]**2.0\
                                            *z[k]**2.0*np.sin(2.0*pi*t[n])

                    term_x[i,j,k,n]     = 2.0*x[i]**2.0*y[j]**2.0*z[k]**2.0\
                                            * np.sin(pi*t[n])\
                                            * np.cos(2.0*pi*t[n])

                    term_y[i,j,k,n]     = 2.0*x[i]**2.0*y[j]**2.0*z[k]**2.0\
                                            * np.cos(pi*t[n])\
                                            * np.cos(2.0*pi*t[n])

                    term_z[i,j,k,n]     = 2.0*x[i]**2.0*y[j]**2.0*z[k]**2.0\
                                            * np.sin(pi*t[n])\
                                            * np.cos(2.0*pi*t[n])
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print("Exact loop --> t_step = %i"      %(n))
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate RHS                                     #
    #---------------------------------------------------------------------#
    RHS_exact       = term_t #+ term_x + term_y + term_z

    print("**** Calculated the exact RHS ****")
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del term_t
    del term_x
    del term_y
    del term_z
    #=====================================================================#
    # Post processing                                                     #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Plotting variables                                                  #
    #---------------------------------------------------------------------#
    zslice      = int((N+1)/2)
    error       = abs(RHS_approx - RHS_exact[:,:,:,t_int])
    print(z[zslice])
    print(t[t_int])
    print(RHS_exact[N,N,zslice,t_int])
    print("Plotting:")
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    [X1, X2]        = np.meshgrid(x,y)
    cnt             = plt.contourf(X1, X2, RHS_approx[:,:,zslice],\
                            500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.title("Approximate RHS")
    plt.colorbar()
    plt.savefig(media_path + "rhs-unit-test-approx.pdf")
    plt.clf()
    print("\tApproximate solution")
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt             = plt.contourf(X1, X2, RHS_exact[:,:,zslice,t_int],\
                            500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.title("Exact RHS")
    plt.colorbar()
    plt.savefig(media_path + "rhs-unit-test-exact.pdf")
    plt.clf()
    print("\tExact solution")
    #---------------------------------------------------------------------#
    # Error solution                                                      #
    #---------------------------------------------------------------------#
    cnt             = plt.contourf(X1, X2, error[:,:,zslice],\
                            500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.title("Error RHS")
    plt.colorbar()
    plt.savefig(media_path + "rhs-unit-test-error.pdf")
    plt.clf()
    print("\tError solution")

    print("**** Successful Run ****")
    sys.exit(0)

#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the P term in the
    kinetic energy transport equation.

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
from strain_rates   import strain_rates
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Terms in the P term                                                     #
#-------------------------------------------------------------------------#
def P_terms(
        U1,                     # velocity-1 component
        U2,                     # velocity-2 component
        U3,                     # velocity-3 component
        Tau11,                  # tau-11 component
        Tau12,                  # tau-12 component
        Tau13,                  # tau-13 component
        Tau22,                  # tau-22 component
        Tau23,                  # tau-23 component
        Tau33,                  # tau-33 component
        h):                     # spatial step size


    """ Calculating the terms in the P term of the kinetic energy transport
    equation """
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33)  = strain_rates(U1, U2, U3, h)
    #---------------------------------------------------------------------#
    # Calculating the terms in the P term                                 #
    #---------------------------------------------------------------------#
    Term1   = np.multiply(Tau11, S11)
    Term2   = 2.0*np.multiply(Tau12, S12)
    Term3   = 2.0*np.multiply(Tau13, S13)
    Term4   = np.multiply(Tau22, S22)
    Term5   = 2.0*np.multiply(Tau23, S23)
    Term6   = np.multiply(Tau33, S33)

    return Term1, Term2, Term3, Term4, Term5, Term6
#-------------------------------------------------------------------------#
# P term                                                                  #
#-------------------------------------------------------------------------#
def p_term(
        U1,                     # velocity-1 component
        U2,                     # velocity-2 component
        U3,                     # velocity-3 component
        Tau11,                  # tau-11 component
        Tau12,                  # tau-12 component
        Tau13,                  # tau-13 component
        Tau22,                  # tau-22 component
        Tau23,                  # tau-23 component
        Tau33,                  # tau-33 component
        h):                     # spatial step size

    """ Calculating the P term """
    #---------------------------------------------------------------------#
    # Calculating the strain rates                                        #
    #---------------------------------------------------------------------#
    (S11, S12, S13, S22, S23, S33)  = strain_rates(U1, U2, U3, h)
    #---------------------------------------------------------------------#
    # Calculating the terms in the P term                                 #
    #---------------------------------------------------------------------#
    Term1   = np.multiply(Tau11, S11)
    Term2   = 2.0*np.multiply(Tau12, S12)
    Term3   = 2.0*np.multiply(Tau13, S13)
    Term4   = np.multiply(Tau22, S22)
    Term5   = 2.0*np.multiply(Tau23, S23)
    Term6   = np.multiply(Tau33, S33)
    #---------------------------------------------------------------------#
    # Calculating the P term                                              #
    #---------------------------------------------------------------------#
    P   = Term1 + Term2 + Term3 + Term4 + Term5 + Term6

    return P
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":
    #=====================================================================#
    # Main preamble                                                       #
    #=====================================================================#
    call(["clear"])
    sep         = os.sep
    pwd         = os.getcwd()
    media_path  = pwd + "%cmedia%c"              %(sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    pi          = np.pi
    N           = 350
    x0          = 0.0
    xf          = 1.0
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x, y)
    dx          = (xf - x0)/N
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating arrays                                                #
    #---------------------------------------------------------------------#
    ux      = np.zeros((N+1, N+1, N+1))
    uy      = np.zeros((N+1, N+1, N+1))
    uz      = np.zeros((N+1, N+1, N+1))
    tau11   = np.zeros((N+1, N+1, N+1))
    tau12   = np.zeros((N+1, N+1, N+1))
    tau13   = np.zeros((N+1, N+1, N+1))
    tau22   = np.zeros((N+1, N+1, N+1))
    tau23   = np.zeros((N+1, N+1, N+1))
    tau33   = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the velocities and the subgrid stresses                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                #---------------------------------------------------------#
                # Velocities                                              #
                #---------------------------------------------------------#
                ux[i,j,k]       = x[i]**2.0*y[j]*z[k]
                uy[i,j,k]       = x[i] + y[j] + z[k]
                uz[i,j,k]       = 3.0*x[i] + 2.0*y[j] + z[k]
                #---------------------------------------------------------#
                # Subgrid stresses                                        #
                #---------------------------------------------------------#
                tau11[i,j,k]    = x[i]**2.0
                tau12[i,j,k]    = x[i]*y[j]
                tau13[i,j,k]    = x[i]*z[k]
                tau22[i,j,k]    = y[j]**2.0
                tau23[i,j,k]    = y[j]*z[k]
                tau33[i,j,k]    = z[k]**2.0
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the approximate solution                                #
    #---------------------------------------------------------------------#
    P_approx    = p_term(ux, uy, uz, tau11, tau12, tau13, tau22, tau23,\
                            tau33, dx)
    #---------------------------------------------------------------------#
    # Deleting variables                                                  #
    #---------------------------------------------------------------------#
    del ux
    del uy
    del uz
    del tau11
    del tau12
    del tau13
    del tau22
    del tau23
    del tau33
    #=====================================================================#
    # Exact solution                                                      #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating arrays                                                #
    #---------------------------------------------------------------------#
    term1   = np.zeros((N+1, N+1, N+1))
    term2   = np.zeros((N+1, N+1, N+1))
    term3   = np.zeros((N+1, N+1, N+1))
    term4   = np.zeros((N+1, N+1, N+1))
    term5   = np.zeros((N+1, N+1, N+1))
    term6   = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating terms in the exact solution                             #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                term1[i,j,k]    = 2.0*x[i]**3.0*y[j]*z[k]
                term2[i,j,k]    = (x[i]*y[j])*(0.5*(x[i]**2.0*z[k]+1))
                term3[i,j,k]    = (x[i]*z[k])*(0.5*(x[i]**2.0*y[j] + 3.0))
                term4[i,j,k]    = y[j]**2.0
                term5[i,j,k]    = (3.0/2.0)*(y[j]*z[k])
                term6[i,j,k]    = z[k]**2.0
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if print_count > 5:
            print(k)
            print_count = 0
        print_count += 1
    #---------------------------------------------------------------------#
    # Calculating the exact solution                                      #
    #---------------------------------------------------------------------#
    P_exact = term1 + 2.0*term2 + 2.0*term3 + term4 + 2.0*term5 + term6
    #=====================================================================#
    # Error                                                               #
    #=====================================================================#
    error   = abs(P_exact - P_approx)
    #=====================================================================#
    # Plotting solutions                                                  #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, P_exact[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "P-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, P_approx[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "P-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Error solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$0 \leq x_{1} \leq 1$")
    plt.ylabel("$0 \leq x_{2} \leq 1$")
    plt.colorbar()
    plt.savefig(media_path + "P-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the C term in the
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
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Terms in the C term                                                     #
#-------------------------------------------------------------------------#
def c_terms(
        U1,                 # velocity-1 component
        U2,                 # velocity-2 component
        U3,                 # velocity-3 component
        Tau11,              # tau-11 component
        Tau12,              # tau-12 component
        Tau13,              # tau-13 component
        Tau22,              # tau-22 component
        Tau23,              # tau-23 component
        Tau33,              # tau-33 component
        h):                 # spatial step size

    """ Calculating the terms in the C term of the kinetic energy transport
    equation """

    #---------------------------------------------------------------------#
    # Calculating the terms in the C term                                 #
    #---------------------------------------------------------------------#
    Term1   = -1.0*np.gradient(Tau11*U1, h, edge_order=2)[0]
    Term2   = -1.0*np.gradient(Tau12*U1, h, edge_order=2)[1]
    Term3   = -1.0*np.gradient(Tau13*U1, h, edge_order=2)[2]
    Term4   = -1.0*np.gradient(Tau12*U2, h, edge_order=2)[0]
    Term5   = -1.0*np.gradient(Tau22*U2, h, edge_order=2)[1]
    Term6   = -1.0*np.gradient(Tau23*U2, h, edge_order=2)[2]
    Term7   = -1.0*np.gradient(Tau13*U3, h, edge_order=2)[0]
    Term8   = -1.0*np.gradient(Tau23*U3, h, edge_order=2)[1]
    Term9   = -1.0*np.gradient(Tau33*U3, h, edge_order=2)[2]

    return Term1, Term2, Term3, Term4, Term5, Term6, Term7, Term8, Term9
#-------------------------------------------------------------------------#
# C term                                                                  #
#-------------------------------------------------------------------------#
def c_term(
        U1,                 # velocity-1 component
        U2,                 # velocity-2 component
        U3,                 # velocity-3 component
        Tau11,              # tau-11 component
        Tau12,              # tau-12 component
        Tau13,              # tau-13 component
        Tau22,              # tau-22 component
        Tau23,              # tau-23 component
        Tau33,              # tau-33 component
        h):                 # spatial step size

    """ Calculating the C term """
    #---------------------------------------------------------------------#
    # Calculating the terms in the C term                                 #
    #---------------------------------------------------------------------#
    Term1   = np.gradient(Tau11*U1, h, edge_order=2)[0]
    Term2   = np.gradient(Tau12*U1, h, edge_order=2)[1]
    Term3   = np.gradient(Tau13*U1, h, edge_order=2)[2]
    Term4   = np.gradient(Tau12*U2, h, edge_order=2)[0]
    Term5   = np.gradient(Tau22*U2, h, edge_order=2)[1]
    Term6   = np.gradient(Tau23*U2, h, edge_order=2)[2]
    Term7   = np.gradient(Tau13*U3, h, edge_order=2)[0]
    Term8   = np.gradient(Tau23*U3, h, edge_order=2)[1]
    Term9   = np.gradient(Tau33*U3, h, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Calculating the C term                                              #
    #---------------------------------------------------------------------#
    C   = -(Term1 + Term2 + Term3 + Term4 + Term5 + Term6 + Term7 + Term8\
                + Term9)
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
    #---------------------------------------------------------------------#
    del Term1
    del Term2
    del Term3
    del Term4
    del Term5
    del Term6
    del Term7
    del Term8
    del Term9

    return C
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
    media_path  = pwd + "%cmedia%c"         %(sep, sep)
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    N           = 400
    x0          = 0.0
    xf          = 1.0
    dx          = (xf-x0)/N
    x           = np.linspace(x0, xf, N+1)
    y           = np.linspace(x0, xf, N+1)
    z           = np.linspace(x0, xf, N+1)
    [X1, X2]    = np.meshgrid(x, y)
    #=====================================================================#
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    ux          = np.zeros((N+1, N+1, N+1))
    uy          = np.zeros((N+1, N+1, N+1))
    uz          = np.zeros((N+1, N+1, N+1))
    tau11       = np.zeros((N+1, N+1, N+1))
    tau12       = np.zeros((N+1, N+1, N+1))
    tau13       = np.zeros((N+1, N+1, N+1))
    tau22       = np.zeros((N+1, N+1, N+1))
    tau23       = np.zeros((N+1, N+1, N+1))
    tau33       = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the velocities and the subgrid stresses                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                ux[i,j,k]   = 2.0*x[i] + y[j] + z[k]
                uy[i,j,k]   = x[i]*y[j]*z[k]
                uz[i,j,k]   = z[k]**2.0
                tau11[i,j,k]    = x[i]**2.0*y[j]*z[k]
                tau12[i,j,k]    = x[i] + y[j] + z[k]
                tau13[i,j,k]    = 2.0*x[i]*z[k]
                tau22[i,j,k]    = y[j]**2.0
                tau23[i,j,k]    = 2.0*x[i]**2.0*y[j]**2.0*z[k]**2.0
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
    C_approx    = c_term(ux, uy, uz, tau11, tau12, tau13, tau22, tau23,\
                            tau33, dx)
    #---------------------------------------------------------------------#
    # Clearing variables                                                  #
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
    # Approximate solution                                                #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    term1   = np.zeros((N+1, N+1, N+1))
    term2   = np.zeros((N+1, N+1, N+1))
    term3   = np.zeros((N+1, N+1, N+1))
    term4   = np.zeros((N+1, N+1, N+1))
    term5   = np.zeros((N+1, N+1, N+1))
    term6   = np.zeros((N+1, N+1, N+1))
    term7   = np.zeros((N+1, N+1, N+1))
    term8   = np.zeros((N+1, N+1, N+1))
    term9   = np.zeros((N+1, N+1, N+1))
    #---------------------------------------------------------------------#
    # Calculating the terms in the exact solution                         #
    #---------------------------------------------------------------------#
    print_count = 0
    for k in range(0, N+1):
        for j in range(0, N+1):
            for i in range(0, N+1):
                #---------------------------------------------------------#
                # terms                                                   #
                #---------------------------------------------------------#
                term1[i,j,k]    = (2.0*x[i] +y[j] +z[k])*(2.0*x[i]*y[j]*z[k]) +\
                                    2.0*(x[i]**2.0*y[j]*z[k])
                term2[i,j,k]    = 3.0*x[i] + 2.0*y[j] + 2.0*z[k]
                term3[i,j,k]    = (2.0*x[i] + y[j] + z[k])*(2.0*x[i]) +\
                                    2.0*x[i]*z[k]
                term4[i,j,k]    = (x[i]*y[j]*z[k]) + (x[i] + y[j] + z[k])*(y[j]*z[k])
                term5[i,j,k]    = 3.0*y[j]**2.0*x[i]*z[k]
                term6[i,j,k]    = 6.0*x[i]**3.0*y[j]**3.0*z[k]**2.0
                term7[i,j,k]    = 2.0*z[k]**3.0
                term8[i,j,k]    = 4.0*x[i]**2.0*y[j]*z[k]**4.0
                term9[i,j,k]    = 4.0*z[k]**3.0
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
    C_exact = -1.0*(term1 + term2 + term3  + term4 + term5 + term6 +\
                        term7 + term8 + term9)
    #=====================================================================#
    # Error                                                               #
    #=====================================================================#
    error   = abs(C_exact - C_approx)
    #=====================================================================#
    # Plotting the solutions                                              #
    #=====================================================================#
    #---------------------------------------------------------------------#
    # Exact solution                                                      #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, C_exact[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"            %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"            %(x0, xf))
    plt.colorbar()
    plt.savefig(media_path + "C-term-exact.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, C_approx[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"            %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"            %(x0, xf))
    plt.colorbar()
    plt.savefig(media_path + "C-term-approx.pdf")
    plt.clf()
    #---------------------------------------------------------------------#
    # Approximate solution                                                #
    #---------------------------------------------------------------------#
    cnt = plt.contourf(X1, X2, error[:,:,int(N/2)], 500, cmap="jet")
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.xlabel("$%.1f \leq x_{1} \leq %.1f$"            %(x0, xf))
    plt.ylabel("$%.1f \leq x_{2} \leq %.1f$"            %(x0, xf))
    plt.colorbar()
    plt.savefig(media_path + "C-term-error.pdf")
    plt.clf()

    print("**** Successful Run ****")
    sys.exit(0)

#!/usr/bin/env python3
"""=============================================================================
Purpose:
    The purpose of this script is to test the user defined functions with the
    results from the simulation for the enstrophy transport equation.

Author:
    Emilio Torres
============================================================================="""
#==============================================================================#
# Preamble                                                                     #
#==============================================================================#
#------------------------------------------------------------------------------#
# Python packages                                                              #
#------------------------------------------------------------------------------#
import os
import sys
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------#
# User packages                                                                #
#------------------------------------------------------------------------------#
from ke_calcs       import ke_average2
#==============================================================================#
# User defined functions                                                       #
#==============================================================================#
#------------------------------------------------------------------------------#
# D/Dt contours                                                                #
#------------------------------------------------------------------------------#
def D_Dt_contours(
        var,
        data,
        norm):

    """ Generating the contours for the (1/Omega)D(Omega)/Dt and (1/K)Dk/Dt """

    data    /= norm
    transport_plots(data, media_path + var + "-" +  str(tval)  , vmin, vmax, dp)

    return
#------------------------------------------------------------------------------#
# Transport term plotting                                                      #
#------------------------------------------------------------------------------#
def countour_generator(
        term,               # term in the transport equation
        var,                # transport variable
        data,               # data being plotted
        norm,               # data to normalize against i.e., Omega/KE
        Zval,               # z-value of interest
        Tval,               # time value of interest
        Vmin,               # minimum plt.colorbar value
        Vmax,               # maximum plt.colorbar value
        DP):                # contour level

    """Subroutine to generate the contour plots for Enstrophy (Omega) and KE"""
    #----------------------------------------------------------------------#
    # Names and paths                                                      #
    #----------------------------------------------------------------------#
    name    = term + '-' + var
    #----------------------------------------------------------------------#
    # Loading A enst data                                                  #
    #----------------------------------------------------------------------#
    data        = data[:,:, zval, tval-tstart]/norm
    transport_plots(data, media_path + name + "-" + str(tval), Vmin, Vmax, DP)
    print(term + '-' + var)

    return
#------------------------------------------------------------------------------#
# Comparison tool                                                              #
#------------------------------------------------------------------------------#
def compare_terms(
        par,                                    # parameter
        sim,                                    # simulation values
        calc,                                   # calculated values
        Valmin,                                 # colorbar minimum value
        Valmax,                                 # colorbar maximum value
        Dp      = False):                       # number of contour values

    """ Tool to compare gradient and spectral fields """
    #--------------------------------------------------------------------------#
    # Default settings                                                         #
    #--------------------------------------------------------------------------#
    if isinstance(par, str) is False:
        par = str(par)
    if Dp is False:
        Dp = (Valmax-Valmin)/500
    #--------------------------------------------------------------------------#
    # Plotting                                                                 #
    #--------------------------------------------------------------------------#
    transport_plots(sim, media_path + par + "-simulation", "Simulation results"\
                    ,Valmin, Valmax, Dp)
    transport_plots(calc, media_path + par + "-calculated", "Calculated results"\
                    ,Valmin, Valmax, Dp)
    print("Plotted " + par)
    #--------------------------------------------------------------------------#
    # Plotting error                                                           #
    #--------------------------------------------------------------------------#
    x           = np.linspace(0,2.0*np.pi,64)
    [X1, X2]    = np.meshgrid(x,x)
    error     = abs(sim - calc)
    cnt         = plt.contourf(X1, X2, error, 500, cmap='jet')
    for c in cnt.collections:
        c.set_edgecolors("face")
    plt.colorbar(cnt)
    plt.xlabel('$0\leq x_{1} \leq 2\pi$')
    plt.ylabel('$0\leq x_{2} \leq 2\pi$')
    plt.title("Error " + par)
    plt.savefig(media_path + par + "-error.png")
    plt.clf()
    print("Error " + par)

    return
#------------------------------------------------------------------------------#
# Plotting tool                                                                #
#------------------------------------------------------------------------------#
def transport_plots(
        var,                    # variable
        plot_name,              # plot name
        Vmin,                   # colorbar minimum
        Vmax,                   # colorbar maximum
        Dp):                    # color resolution

    """ Subroutine to plot transport terms """
    #--------------------------------------------------------------------------#
    # Domain variables                                                         #
    #--------------------------------------------------------------------------#
    pi          = np.pi
    x           = np.linspace(0,2.0*pi,64)
    [X1, X2]    = np.meshgrid(x,x)
    #--------------------------------------------------------------------------#
    # Plotting                                                                 #
    #--------------------------------------------------------------------------#
    fig0, ax0   = plt.subplots(1,1,)
    cnt         = ax0.contourf(X1, X2, var, np.arange(Vmin, Vmax, Dp),\
                            cmap='jet', extend='both')
    for c in cnt.collections:
        c.set_edgecolors("face")
    cbar0   = plt.colorbar(cnt)
    plt.xlabel('$0\leq x_{1} \leq 2\pi$')
    plt.ylabel('$0\leq x_{2} \leq 2\pi$')
    plt.savefig(plot_name + ".png")
    plt.clf()

    return
#==============================================================================#
# Main                                                                         #
#==============================================================================#
if __name__ == "__main__":
    #--------------------------------------------------------------------------#
    # Main preamble                                                            #
    #--------------------------------------------------------------------------#
    call(["clear"])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + "%c..%cdata%c"              %(sep, sep, sep)
    media_path  = pwd + "%c..%cmedia%c"             %(sep, sep, sep)
    #--------------------------------------------------------------------------#
    # Domain variables                                                         #
    #--------------------------------------------------------------------------#
    zval    = 31
    tstart  = 1337
    times   = [tstart, 1357, 1377, 1397, 1417, 1437, 1457, 1477, 1497, 1517,\
                1537, 1557, 1577, 1597, 1617, 1637, 1657, 1677, 1697, 1717,\
                1737, 1757, 1777, 1797, 1817, 1837, 1857, 1877, 1897, 1917,\
                1937]
    #--------------------------------------------------------------------------#
    # Loading field data                                                       #
    #--------------------------------------------------------------------------#
    print("Loading field data")
    keavg   = np.loadtxt("../ke-avg.txt", skiprows=2)
    time    = keavg[:,1]
    keavg   = keavg[:,2]
    print("\ttime")
    print("\tke-average")
    ke      = np.load(data_path + "ke.npy")
    print("\tkinetic energy")
    enst    = np.load(data_path + "enst.npy")
    print("\tenstrophy")
    #--------------------------------------------------------------------------#
    # Loading data                                                             #
    #--------------------------------------------------------------------------#
    print("Loading enstrophy data")
    Aenst   = np.load(data_path + "A-" + "enst" + ".npy")
    print("\tA")
    Benst   = np.load(data_path + "B-" + "enst" + ".npy")
    print("\tB")
    Denst   = np.load(data_path + "D-" + "enst" + ".npy")
    print("\tD")
    Pienst  = np.load(data_path + "Pi-" + "enst" + ".npy")
    print("\tPi")
    Penst   = np.load(data_path + "P-" + "enst" + ".npy")
    print("\tP")
    #--------------------------------------------------------------------------#
    # Loading data                                                             #
    #--------------------------------------------------------------------------#
    print("Loading kinetic energy data")
    Ake     = np.load(data_path + "A-" + "ke" + ".npy")
    print("\tA")
    Bke     = np.load(data_path + "B-" + "ke" + ".npy")
    print("\tB")
    Dke     = np.load(data_path + "D-" + "ke" + ".npy")
    print("\tD")
    Cke     = np.load(data_path + "C-" + "ke" + ".npy")
    print("\tC")
    Pke     = np.load(data_path + "P-" + "ke" + ".npy")
    print("\tP")
    print("Loaded KE data")
    #--------------------------------------------------------------------------#
    # Enstrophy flags                                                          #
    #--------------------------------------------------------------------------#
    enstavg_flag    = False
    keavg_flag      = False
    a_enst_flag     = False
    b_enst_flag     = False
    d_enst_flag     = False
    pi_enst_flag    = False
    p_enst_flag     = False
    enst_flag       = False
    enst2_flag      = False
    #--------------------------------------------------------------------------#
    # kinetic energy flags                                                     #
    #--------------------------------------------------------------------------#
    a_ke_flag       = True
    b_ke_flag       = True
    d_ke_flag       = True
    c_ke_flag       = True
    p_ke_flag       = True
    ke_flag         = True
    ke2_flag        = True
    #==========================================================================#
    # KE average                                                               #
    #==========================================================================#
    for count, tval in enumerate(times):
        print(time[tval+1])
        print(count)
        print(tval)
        if enstavg_flag is True:
            enstavg   = np.load(data_path + "enst-avg.npy")
            plt.plot(time[0:tval+1], enstavg[0:tval+1], 'r', lw=1.5)
            plt.xlabel("Simulation time")
            plt.ylabel("$\\langle \Omega \\rangle$")
            plt.grid(True)
            plt.savefig(media_path + "enst-average" + str(tval) + ".png")
            plt.clf()
            print("Plotted Omega average")
        if keavg_flag is True:
            plt.plot(time[0:tval+1], keavg[0:tval+1], 'r', lw=1.5)
            plt.xlabel("Simulation time")
            plt.ylabel("$\\langle k \\rangle$")
            plt.grid(True)
            plt.savefig(media_path + "ke-average" + str(tval) + ".png")
            plt.clf()
            print("Plotted KE average")
        #==========================================================================#
        # Enstrophy transport                                                      #
        #==========================================================================#
        #----------------------------------------------------------------------#
        # Plotting settings                                                    #
        #----------------------------------------------------------------------#
        c1          = 80
        vmin        = -c1
        vmax        = c1
        dp          = (vmax - vmin)/500
        enst_norm   = enst[:,:,zval,tval-tstart]
        if a_enst_flag is True:
            countour_generator("A","enst", Aenst, enst_norm, zval, tval-tstart, vmin, vmax, dp)
        if b_enst_flag is True:
            countour_generator("B","enst", Benst, enst_norm, zval, tval-tstart, vmin, vmax, dp)
        if d_enst_flag is True:
            countour_generator("D","enst", Denst, enst_norm, zval, tval-tstart, vmin, vmax, dp)
        if pi_enst_flag is True:
            countour_generator("Pi","enst", Pienst, enst_norm, zval, tval-tstart, vmin, vmax, dp)
        if p_enst_flag is True:
            countour_generator("P","enst", Penst, enst_norm, zval, tval-tstart, vmin, vmax, dp)
        if enst_flag is True:
            A       = Aenst[:,:, zval, tval-tstart]
            B       = Benst[:,:, zval, tval-tstart]
            D       = Denst[:,:, zval, tval-tstart]
            trans   = Pienst[:,:, zval, tval-tstart]
            P       = Penst[:,:, zval, tval-tstart]
            D_Dt_contours("enst", A + B + D + trans + P, enst_norm)
        if enst2_flag is True:
            vmin        = 2.0*np.std(enst[:,:,zval,tval-tstart])
            vmax        = np.amax(enst[:,:,zval,tval-tstart])
            dp          = (vmax-vmin)/500
            D_Dt_contours("enst-2", enst[:,:,zval,tval-tstart], np.ones((64,64)))
        #==========================================================================#
        # Kinetic energy transport                                                 #
        #==========================================================================#
        #----------------------------------------------------------------------#
        # Plotting settings                                                    #
        #----------------------------------------------------------------------#
        c1          = 25
        vmin        = -c1
        vmax        = c1
        dp          = (vmax - vmin)/500
        ke_norm     = ke[:,:,zval,tval-tstart]
        if a_ke_flag is True:
            countour_generator("A","ke", Ake, ke_norm, zval, tval-tstart, vmin, vmax, dp)
        if b_ke_flag is True:
            countour_generator("B","ke", Bke, ke_norm, zval, tval-tstart, vmin, vmax, dp)
        if d_ke_flag is True:
            countour_generator("D","ke", Dke, ke_norm, zval, tval-tstart, vmin, vmax, dp)
        if c_ke_flag is True:
            countour_generator("C","ke", Cke, ke_norm, zval, tval-tstart, vmin, vmax, dp)
        if p_ke_flag is True:
            countour_generator("P","ke", Pke, ke_norm, zval, tval-tstart, vmin, vmax, dp)
        if ke_flag is True:
            c1          = 25
            vmin        = -c1
            vmax        = c1
            dp          = (vmax-vmin)/500
            A       = Ake[:,:, zval, tval-tstart]
            B       = Bke[:,:, zval, tval-tstart]
            D       = Dke[:,:, zval, tval-tstart]
            trans   = Cke[:,:, zval, tval-tstart]
            P       = Pke[:,:, zval, tval-tstart]
            D_Dt_contours("ke", A + B + D + trans + P, ke_norm)
        if ke2_flag is True:
            vmin    = 2.0*np.std(ke[:,:,zval,tval-tstart]) 
            vmax    = np.amax(ke[:,:,zval,tval-tstart]) 
            dp      = (vmax-vmin)/500
            D_Dt_contours("ke-2", ke[:,:,zval,tval-tstart], np.ones((64,64)))

    print("**** successful run ****")
    sys.exit(0)

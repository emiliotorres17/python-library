#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is extract the data for BS run.

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
#-------------------------------------------------------------------------#
# User defined packages                                                   #
#-------------------------------------------------------------------------#
from npy_extract    import npy_time
from npy_extract    import npy_velocity_general
from npy_extract    import npy_vorticity_general
from npy_extract    import npy_transport_term_general
from npy_extract    import npy_S_general
from npy_extract    import npy_tau_general
from npy_extract    import npy_psi_general
from ke_calcs       import ke_average2
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
    data_path       = pwd + "%c..%cdata%c"              %(sep, sep, sep)
    time_path       = pwd + "%c..%ctime%c"              %(sep, sep, sep)
    A_path_enst     = pwd + "%c..%cA-enst%c"            %(sep, sep, sep)
    B_path_enst     = pwd + "%c..%cB-enst%c"            %(sep, sep, sep)
    D_path_enst     = pwd + "%c..%cD-enst%c"            %(sep, sep, sep)
    C_path_enst     = pwd + "%c..%cC-enst%c"            %(sep, sep, sep)
    P_path_enst     = pwd + "%c..%cP-enst%c"            %(sep, sep, sep)
    A_path_ke       = pwd + "%c..%cA-ke%c"              %(sep, sep, sep)
    B_path_ke       = pwd + "%c..%cB-ke%c"              %(sep, sep, sep)
    D_path_ke       = pwd + "%c..%cD-ke%c"              %(sep, sep, sep)
    C_path_ke       = pwd + "%c..%cC-ke%c"              %(sep, sep, sep)
    P_path_ke       = pwd + "%c..%cP-ke%c"              %(sep, sep, sep)
    ke_path         = pwd + "%c..%cke%c"                %(sep, sep, sep)
    enst_path       = pwd + "%c..%censt%c"              %(sep, sep, sep)
    vel1_path       = pwd + "%c..%cvelocity1%c"         %(sep, sep, sep)
    vel2_path       = pwd + "%c..%cvelocity2%c"         %(sep, sep, sep)
    vel3_path       = pwd + "%c..%cvelocity3%c"         %(sep, sep, sep)
    omega1_path     = pwd + "%c..%comega1%c"            %(sep, sep, sep)
    omega2_path     = pwd + "%c..%comega2%c"            %(sep, sep, sep)
    omega3_path     = pwd + "%c..%comega3%c"            %(sep, sep, sep)
    S_path          = pwd + "%c..%cstrain-rates%c"      %(sep, sep, sep)
    tau_path        = pwd + "%c..%ctau%c"               %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    tstart      = 1337
    tf          = 2337
    time_flag   = False
    enst_flag   = False
    ke_flag     = False
    fields_flag = True
    avg_flag    = True
    vel_flag    = False
    omega_flag  = False
    #=====================================================================#
    # Time                                                                #
    #=====================================================================#
    if time_flag is True:
        #-----------------------------------------------------------------#
        # Extracting time                                                 #
        #-----------------------------------------------------------------#
        time        = npy_time(tf, time_path)
        np.save(data_path + "time.npy", time)
        del time
    #=====================================================================#
    # Enstrophy transport terms                                           #
    #=====================================================================#
    if enst_flag is True:
        #-----------------------------------------------------------------#
        # Extracting A spectral                                           #
        #-----------------------------------------------------------------#
        A   = npy_transport_term_general('A', 64, 16, tstart, tf, A_path_enst)
        np.save(data_path + "A-enst.npy", A)
        del A
        #-----------------------------------------------------------------#
        # Extracting B spectral                                           #
        #-----------------------------------------------------------------#
        B   = npy_transport_term_general('B', 64, 16, tstart, tf, B_path_enst)
        np.save(data_path + "B-enst.npy", B)
        del B
        #-----------------------------------------------------------------#
        # Extracting D spectral                                           #
        #-----------------------------------------------------------------#
        D   = npy_transport_term_general('D', 64, 16, tstart, tf, D_path_enst)
        np.save(data_path + "D-enst.npy", D)
        del D
        #-----------------------------------------------------------------#
        # Extracting Pi spectral                                          #
        #-----------------------------------------------------------------#
        C   = npy_transport_term_general('C', 64, 16, tstart, tf, C_path_enst)
        np.save(data_path + "C-enst.npy", D)
        del C
        #-----------------------------------------------------------------#
        # Extracting P spectral                                           #
        #-----------------------------------------------------------------#
        P   = npy_transport_term_general('P', 64, 16, tstart, tf, P_path_enst)
        np.save(data_path + "P-enst.npy", P)
        del P
    #=====================================================================#
    # KE transport terms                                                  #
    #=====================================================================#
    if ke_flag is True:
        #-----------------------------------------------------------------#
        # Extracting A spectral                                           #
        #-----------------------------------------------------------------#
        A   = npy_transport_term_general('A', 64, 16, tstart, tf, A_path_ke)
        np.save(data_path + "A-ke.npy", A)
        del A
        #-----------------------------------------------------------------#
        # Extracting B spectral                                           #
        #-----------------------------------------------------------------#
        B   = npy_transport_term_general('B', 64, 16, tstart, tf, B_path_ke)
        np.save(data_path + "B-ke.npy", B)
        del B
        #-----------------------------------------------------------------#
        # Extracting D spectral                                           #
        #-----------------------------------------------------------------#
        D   = npy_transport_term_general('D', 64, 16, tstart, tf, D_path_ke)
        np.save(data_path + "D-ke.npy", D)
        del D
        #-----------------------------------------------------------------#
        # Extracting Pi spectral                                          #
        #-----------------------------------------------------------------#
        C   = npy_transport_term_general('C', 64, 16, tstart, tf, C_path_ke)
        np.save(data_path + "C-ke.npy", C)
        del C
        #-----------------------------------------------------------------#
        # Extracting P spectral                                           #
        #-----------------------------------------------------------------#
        P   = npy_transport_term_general('P', 64, 16, tstart, tf, P_path_ke)
        np.save(data_path + "P-ke.npy", P)
        del P
    #=====================================================================#
    # Enstrophy and kinetic energy fields                                 #
    #=====================================================================#
    if fields_flag is True:
        #-----------------------------------------------------------------#
        # KE fields                                                       #
        #-----------------------------------------------------------------#
        ke      = npy_transport_term_general('ke', 64, 16, tstart, tf, ke_path)
        np.save(data_path + "ke.npy", ke)
        del ke
        #-----------------------------------------------------------------#
        # Enstrophy fields                                                #
        #-----------------------------------------------------------------#
        enst    = npy_transport_term_general('enstrophy', 64, 16, tstart, tf, enst_path)
        np.save(data_path + "enst.npy", enst)
        del enst
    #=====================================================================#
    # Velocity fields                                                     #
    #=====================================================================#
    if vel_flag is True:
        #-----------------------------------------------------------------#
        # velocity-1                                                      #
        #-----------------------------------------------------------------#
        vel1    = npy_velocity_general(1, 64, 16, 0 tf, vel1_path)
        np.save(data_path + 'velocity1.npy', vel1)
        del vel1
        #-----------------------------------------------------------------#
        # velocity-2                                                      #
        #-----------------------------------------------------------------#
        vel2    = npy_velocity_general(2, 64, 16, 0 tf, vel2_path)
        np.save(data_path + 'velocity2.npy', vel2)
        del vel2
        #-----------------------------------------------------------------#
        # velocity-3                                                      #
        #-----------------------------------------------------------------#
        vel3    = npy_velocity_general(3, 64, 16, 0 tf, vel3_path)
        np.save(data_path + 'velocity3.npy', vel3)
        del vel3
    #=====================================================================#
    # Vorticity fields                                                    #
    #=====================================================================#
    if omega_flag is True:
        #-----------------------------------------------------------------#
        # velocity-1                                                      #
        #-----------------------------------------------------------------#
        omega1  = npy_vorticity_general(1, 64, 16, 0 tf, omega1_path)
        np.save(data_path + 'omega1.npy', omega1)
        del omega1
        #-----------------------------------------------------------------#
        # velocity-2                                                      #
        #-----------------------------------------------------------------#
        omega2  = npy_vorticity_general(2, 64, 16, 0 tf, omega2_path)
        np.save(data_path + 'omega2.npy', omega2)
        del omega2
        #-----------------------------------------------------------------#
        # velocity-3                                                      #
        #-----------------------------------------------------------------#
        omega3  = npy_vorticity_general(3, 64, 16, 0 tf, omega3_path)
        np.save(data_path + 'omega3.npy', omega3)
        del omega3
    #=====================================================================#
    # Strain rates                                                        #
    #=====================================================================#
    if S_flag is True:
        S   = np.zeros((6,64,64,64,tf+1))
        #-----------------------------------------------------------------#
        # S-11                                                            #
        #-----------------------------------------------------------------#
        s11     = npy_S_general(11, 64, 16, 0, tf, S_path)
        S[0]    = s11
        #np.save(data_path + 's11.npy', s11)
        del s11
        #-----------------------------------------------------------------#
        # S-12                                                            #
        #-----------------------------------------------------------------#
        s12     = npy_S_general(12, 64, 16, 0, tf, S_path)
        S[1]    = s12
        #np.save(data_path + 's12.npy', s12)
        del s12
        #-----------------------------------------------------------------#
        # S-13                                                            #
        #-----------------------------------------------------------------#
        s13     = npy_S_general(13, 64, 16, 0, tf, S_path)
        S[2]    = s13
        #np.save(data_path + 's13.npy', s13)
        del s13
        #-----------------------------------------------------------------#
        # S-22                                                            #
        #-----------------------------------------------------------------#
        s22     = npy_S_general(22, 64, 16, 0, tf, S_path)
        S[3]    = s22
        #np.save(data_path + 's22.npy', s22)
        del s22
        #-----------------------------------------------------------------#
        # S-23                                                            #
        #-----------------------------------------------------------------#
        s23     = npy_S_general(23, 64, 16, 0, tf, S_path)
        S[4]    = s23
        #np.save(data_path + 's23.npy', s23)
        del s23
        #-----------------------------------------------------------------#
        # S-33                                                            #
        #-----------------------------------------------------------------#
        s33 = npy_S_general(33, 64, 16, 0, tf, S_path)
        S[5]    = s23
        #np.save(data_path + 's33.npy', s33)
        del s33
        np.save(data_path + 'strain-rates.npy', S)
    #=====================================================================#
    # Subgrid stress                                                      #
    #=====================================================================#
    if tau_flag is True:
        tau   = np.zeros((6,64,64,64,tf+1))
        #-----------------------------------------------------------------#
        # tau-11                                                          #
        #-----------------------------------------------------------------#
        tau11     = npy_tau_general(11, 64, 16, 0, tf, S_path)
        tau[0]    = tau11
        del tau11
        #-----------------------------------------------------------------#
        # tau-12                                                          #
        #-----------------------------------------------------------------#
        tau12     = npy_tau_general(12, 64, 16, 0, tf, S_path)
        tau[1]    = tau12
        del tau12
        #-----------------------------------------------------------------#
        # tau-13                                                          #
        #-----------------------------------------------------------------#
        tau13     = npy_tau_general(13, 64, 16, 0, tf, S_path)
        tau[2]    = tau13
        del tau13
        #-----------------------------------------------------------------#
        # tau-22                                                          #
        #-----------------------------------------------------------------#
        tau22     = npy_tau_general(22, 64, 16, 0, tf, S_path)
        tau[3]    = tau22
        del tau22
        #-----------------------------------------------------------------#
        # tau-23                                                          #
        #-----------------------------------------------------------------#
        tau23     = npy_tau_general(23, 64, 16, 0, tf, S_path)
        tau[4]    = tau23
        del tau23
        #-----------------------------------------------------------------#
        # tau-33                                                          #
        #-----------------------------------------------------------------#
        tau33   = npy_tau_general(33, 64, 16, 0, tf, S_path)
        tau[5]  = tau23
        del tau33
        np.save(data_path + 'tau.npy', tau)

    print("**** Successful Run ****")
    sys.exit(0)

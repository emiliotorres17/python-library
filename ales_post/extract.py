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
from npy_extract    import npy_time_interval
from npy_extract    import npy_velocity_general
from npy_extract    import npy_vorticity_general
from npy_extract    import npy_transport_term_general
from npy_extract    import npy_S_general
from npy_extract    import npy_tau_general
from npy_extract    import npy_psi_general
from ke_calcs       import ke_average2
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Deleting variables                                                      #
#-------------------------------------------------------------------------#
def delete_var(
        path,
        name):

    """ Deleting variables """
    #---------------------------------------------------------------------#
    # Navigating to path                                                  #
    #---------------------------------------------------------------------#
    os.chdir(path)
    #---------------------------------------------------------------------#
    # Looping over files                                                  #
    #---------------------------------------------------------------------#
    count = 0
    for i, filename in enumerate(os.listdir()):
        if filename.endswith('.npy'):
           os.remove(filename)
        #-----------------------------------------------------------------#
        # Print statement                                                 #
        #-----------------------------------------------------------------#
        if count > 100:
            print('Deleting ' + name + ' time step --> %i'       %(i)) 
            count = 0
        count += 1
    os.chdir(pwd)

    return
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
    omega1_path     = pwd + "%c..%cvorticity1%c"        %(sep, sep, sep)
    omega2_path     = pwd + "%c..%cvorticity2%c"        %(sep, sep, sep)
    omega3_path     = pwd + "%c..%cvorticity3%c"        %(sep, sep, sep)
    tau_path        = pwd + "%c..%ctau%c"               %(sep, sep, sep)
    S_path          = pwd + "%c..%cstrain-rates%c"      %(sep, sep, sep)
    cs2_path        = pwd + "%c..%ccs2%c"               %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    tstart      = 0
    tf          = 503
    num_proc    = 16
    N           = 64
    del_flag    = False
    time_flag   = True
    enst_flag   = True
    ke_flag     = True
    fields_flag = True
    avg_flag    = False
    vel_flag    = False
    omega_flag  = False
    S_flag      = False
    tau_flag    = False
    cs2_flag    = False
    #=====================================================================#
    # Time                                                                #
    #=====================================================================#
    if time_flag is True:
        #-----------------------------------------------------------------#
        # Extracting time                                                 #
        #-----------------------------------------------------------------#
        time        = npy_time_interval(tstart, tf, time_path)
        np.save(data_path + "time.npy", time)
        if del_flag is True:
            delete_var(time_path, 'time')
        del time 
    #=====================================================================#
    # Enstrophy transport terms                                           #
    #=====================================================================#
    if enst_flag is True:
        #-----------------------------------------------------------------#
        # Extracting A spectral                                           #
        #-----------------------------------------------------------------#
        A   = npy_transport_term_general('A', N, num_proc, tstart, tf, A_path_enst)
        np.save(data_path + "A-enst.npy", A)
        if del_flag is True:
            delete_var(A_path_enst, 'A-enst')
        del A
        #-----------------------------------------------------------------#
        # Extracting B spectral                                           #
        #-----------------------------------------------------------------#
        B   = npy_transport_term_general('B', N, num_proc, tstart, tf, B_path_enst)
        np.save(data_path + "B-enst.npy", B)
        if del_flag is True:
            delete_var(B_path_enst, 'B-enst')
        del B
        #-----------------------------------------------------------------#
        # Extracting D spectral                                           #
        #-----------------------------------------------------------------#
        D   = npy_transport_term_general('D', N, num_proc, tstart, tf, D_path_enst)
        np.save(data_path + "D-enst.npy", D)
        if del_flag is True:
            delete_var(D_path_enst, 'D-enst')
        del D
        #-----------------------------------------------------------------#
        # Extracting C spectral                                           #
        #-----------------------------------------------------------------#
        C   = npy_transport_term_general('C', N, num_proc, tstart, tf, C_path_enst)
        np.save(data_path + "C-enst.npy", C)
        if del_flag is True:
            delete_var(C_path_enst, 'C-enst')
        del C
        #-----------------------------------------------------------------#
        # Extracting P spectral                                           #
        #-----------------------------------------------------------------#
        P   = npy_transport_term_general('P', N, num_proc, tstart, tf, P_path_enst)
        np.save(data_path + "P-enst.npy", P)
        if del_flag is True:
            delete_var(P_path_enst, 'P-enst')
        del P
    #=====================================================================#
    # KE transport terms                                                  #
    #=====================================================================#
    if ke_flag is True:
        #-----------------------------------------------------------------#
        # Extracting A spectral                                           #
        #-----------------------------------------------------------------#
        A   = npy_transport_term_general('A', N, num_proc, tstart, tf, A_path_ke)
        np.save(data_path + "A-ke.npy", A)
        if del_flag is True:
            delete_var(A_path_ke, 'A-ke')
        del A
        #-----------------------------------------------------------------#
        # Extracting B spectral                                           #
        #-----------------------------------------------------------------#
        B   = npy_transport_term_general('B', N, num_proc, tstart, tf, B_path_ke)
        np.save(data_path + "B-ke.npy", B)
        if del_flag is True:
            delete_var(B_path_ke, 'B-ke')
        del B
        #-----------------------------------------------------------------#
        # Extracting D spectral                                           #
        #-----------------------------------------------------------------#
        D   = npy_transport_term_general('D', N, num_proc, tstart, tf, D_path_ke)
        np.save(data_path + "D-ke.npy", D)
        if del_flag is True:
            delete_var(D_path_ke, 'D-ke')
        del D
        #-----------------------------------------------------------------#
        # Extracting Pi spectral                                          #
        #-----------------------------------------------------------------#
        C   = npy_transport_term_general('C', N, num_proc, tstart, tf, C_path_ke)
        np.save(data_path + "C-ke.npy", C)
        if del_flag is True:
            delete_var(C_path_ke, 'C-ke')
        del C
        #-----------------------------------------------------------------#
        # Extracting P spectral                                           #
        #-----------------------------------------------------------------#
        P   = npy_transport_term_general('P', N, num_proc, tstart, tf, P_path_ke)
        np.save(data_path + "P-ke.npy", P)
        if del_flag is True:
            delete_var(P_path_ke, 'P-ke')
        del P
    #=====================================================================#
    # Enstrophy and kinetic energy fields                                 #
    #=====================================================================#
    if fields_flag is True:
        #-----------------------------------------------------------------#
        # KE fields                                                       #
        #-----------------------------------------------------------------#
        ke      = npy_transport_term_general('ke', N, num_proc, tstart, tf, ke_path)
        np.save(data_path + "ke.npy", ke)
        if del_flag is True:
            delete_var(ke_path, 'ke')
        del ke
        #-----------------------------------------------------------------#
        # Enstrophy fields                                                #
        #-----------------------------------------------------------------#
        enst    = npy_transport_term_general('enstrophy', N, num_proc, tstart, tf, enst_path)
        np.save(data_path + "enst.npy", enst)
        if del_flag is True:
            delete_var(enst_path, 'enst')
        del enst
    #=====================================================================#
    # Velocity fields                                                     #
    #=====================================================================#
    if vel_flag is True:
        #-----------------------------------------------------------------#
        # velocity-1                                                      #
        #-----------------------------------------------------------------#
        vel1    = npy_velocity_general(1, N, num_proc, tstart, tf, vel1_path)
        np.save(data_path + 'velocity1.npy', vel1)
        if del_flag is True:
            delete_var(vel1_path, 'vel-1')
        del vel1
        #-----------------------------------------------------------------#
        # velocity-2                                                      #
        #-----------------------------------------------------------------#
        vel2    = npy_velocity_general(2, N, num_proc, tstart, tf, vel2_path)
        np.save(data_path + 'velocity2.npy', vel2)
        if del_flag is True:
            delete_var(vel2_path, 'vel-2')
        del vel2
        #-----------------------------------------------------------------#
        # velocity-3                                                      #
        #-----------------------------------------------------------------#
        vel3    = npy_velocity_general(3, N, num_proc, tstart, tf, vel3_path)
        np.save(data_path + 'velocity3.npy', vel3)
        if del_flag is True:
            delete_var(vel3_path, 'vel-3')
        del vel3
    #=====================================================================#
    # Vorticity fields                                                    #
    #=====================================================================#
    if omega_flag is True:
        #-----------------------------------------------------------------#
        # velocity-1                                                      #
        #-----------------------------------------------------------------#
        omega1  = npy_vorticity_general(1, N, num_proc, tstart, tf, omega1_path)
        np.save(data_path + 'omega1.npy', omega1)
        if del_flag is True:
            delete_var(omega1_path, 'omega-1')
        del omega1
        #-----------------------------------------------------------------#
        # velocity-2                                                      #
        #-----------------------------------------------------------------#
        omega2  = npy_vorticity_general(2, N, num_proc, tstart, tf, omega2_path)
        np.save(data_path + 'omega2.npy', omega2)
        if del_flag is True:
            delete_var(omega2_path, 'omega-2')
        del omega2
        #-----------------------------------------------------------------#
        # velocity-3                                                      #
        #-----------------------------------------------------------------#
        omega3  = npy_vorticity_general(3, N, num_proc, tstart, tf, omega3_path)
        np.save(data_path + 'omega3.npy', omega3)
        if del_flag is True:
            delete_var(omega3_path, 'omega-3')
        del omega3
    #=====================================================================#
    # Strain rates                                                        #
    #=====================================================================#
    if S_flag is True:
        #-----------------------------------------------------------------#
        # S-11                                                            #
        #-----------------------------------------------------------------#
        s11     = npy_S_general(11, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's11.npy', s11)
        del s11
        #-----------------------------------------------------------------#
        # S-12                                                            #
        #-----------------------------------------------------------------#
        s12 = npy_S_general(12, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's12.npy', s12)
        del s12
        #-----------------------------------------------------------------#
        # S-13                                                            #
        #-----------------------------------------------------------------#
        s13 = npy_S_general(13, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's13.npy', s13)
        del s13
        #-----------------------------------------------------------------#
        # S-22                                                            #
        #-----------------------------------------------------------------#
        s22 = npy_S_general(22, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's22.npy', s22)
        del s22
        #-----------------------------------------------------------------#
        # S-23                                                            #
        #-----------------------------------------------------------------#
        s23 = npy_S_general(23, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's23.npy', s23)
        del s23
        #-----------------------------------------------------------------#
        # S-33                                                            #
        #-----------------------------------------------------------------#
        s33 = npy_S_general(33, N, num_proc, tstart, tf, S_path)
        np.save(data_path + 's33.npy', s33)
        del s33
        if del_flag is True:
            delete_var(S_path, 'strain-rates')
    #=====================================================================#
    # Subgrid stress                                                      #
    #=====================================================================#
    if tau_flag is True:
        #-----------------------------------------------------------------#
        # tau-11                                                          #
        #-----------------------------------------------------------------#
        tau11   = npy_tau_general(11, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau11.npy', tau11)
        del tau11
        #-----------------------------------------------------------------#
        # tau-12                                                          #
        #-----------------------------------------------------------------#
        tau12   = npy_tau_general(12, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau12.npy', tau12)
        del tau12
        #-----------------------------------------------------------------#
        # tau-13                                                          #
        #-----------------------------------------------------------------#
        tau13   = npy_tau_general(13, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau13.npy', tau13)
        del tau13
        #-----------------------------------------------------------------#
        # tau-22                                                          #
        #-----------------------------------------------------------------#
        tau22   = npy_tau_general(22, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau22.npy', tau22)
        del tau22
        #-----------------------------------------------------------------#
        # tau-23                                                          #
        #-----------------------------------------------------------------#
        tau23   = npy_tau_general(23, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau23.npy', tau23)
        del tau23
        #-----------------------------------------------------------------#
        # tau-33                                                          #
        #-----------------------------------------------------------------#
        tau33   = npy_tau_general(33, N, num_proc, tstart, tf, tau_path)
        np.save(data_path + 'tau33.npy', tau33)
        del tau33
        if del_flag is True:
            delete_var(tau_path, 'subgrid-stress')
    #=====================================================================#
    # Dynamic Smagorinsky                                                 #
    #=====================================================================#
    if cs2_flag is True:
        cs2 = npy_transport_term_general('Cs2', N, num_proc, tstart, tf, cs2_path)
        np.save(data_path + 'cs2.npy', cs2)
        del cs2
        if del_flag is True:
            delete_var(cs2_path, 'cs2')
    print("**** Successful Run ****")
    sys.exit(0)

#!/usr/bin/env python3
"""=============================================================================
Purpose:
    The purpose of this script is extract the data for BS run.

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
#------------------------------------------------------------------------------#
# User defined packages                                                        #
#------------------------------------------------------------------------------#
from npy_extract    import npy_time
from npy_extract    import npy_velocity_general
from npy_extract    import npy_vorticity_general
from npy_extract    import npy_transport_term_general
from npy_extract    import npy_S_general
from npy_extract    import npy_tau_general
from npy_extract    import npy_psi_general
from ke_calcs       import ke_average2
#==============================================================================#
# Main                                                                         #
#==============================================================================#
if __name__ == "__main__":
    #--------------------------------------------------------------------------#
    # Main preamble                                                            #
    #--------------------------------------------------------------------------#
    call(["clear"])
    sep             = os.sep
    pwd             = os.getcwd()
    data_path       = pwd + "%c..%cdata%c"              %(sep, sep, sep)
    time_path       = pwd + "%c..%ctime%c"              %(sep, sep, sep)
    A_path_enst     = pwd + "%c..%cA-enst%c"            %(sep, sep, sep)
    B_path_enst     = pwd + "%c..%cB-enst%c"            %(sep, sep, sep)
    D_path_enst     = pwd + "%c..%cD-enst%c"            %(sep, sep, sep)
    Pi_path_enst    = pwd + "%c..%cPi-enst%c"           %(sep, sep, sep)
    P_path_enst     = pwd + "%c..%cP-enst%c"            %(sep, sep, sep)
    A_path_ke       = pwd + "%c..%cA-ke%c"              %(sep, sep, sep)
    B_path_ke       = pwd + "%c..%cB-ke%c"              %(sep, sep, sep)
    D_path_ke       = pwd + "%c..%cD-ke%c"              %(sep, sep, sep)
    C_path_ke       = pwd + "%c..%cC-ke%c"              %(sep, sep, sep)
    P_path_ke       = pwd + "%c..%cP-ke%c"              %(sep, sep, sep)
    ke_path         = pwd + "%c..%cke%c"                %(sep, sep, sep)
    enst_path       = pwd + "%c..%censt%c"              %(sep, sep, sep)
    #--------------------------------------------------------------------------#
    # Domain variables                                                         #
    #--------------------------------------------------------------------------#
    tstart      = 1337
    tf          = 2337
    time_flag   = False
    enst_flag   = False
    ke_flag     = False
    fields_flag = True
    avg_flag    = True
    #==========================================================================#
    # Time                                                                     #
    #==========================================================================#
    if time_flag is True:
        #----------------------------------------------------------------------#
        # Extracting time                                                      #
        #----------------------------------------------------------------------#
        time        = npy_time(tf, time_path)
        np.save(data_path + "time.npy", time)
        del time
    #==========================================================================#
    # Enstrophy transport terms                                                #
    #==========================================================================#
    if enst_flag is True:
        #----------------------------------------------------------------------#
        # Extracting A spectral                                                #
        #----------------------------------------------------------------------#
        A   = npy_transport_term_general('A', 64, 16, tstart, tf, A_path_enst)
        np.save(data_path + "A-enst.npy", A)
        del A
        #----------------------------------------------------------------------#
        # Extracting B spectral                                                #
        #----------------------------------------------------------------------#
        B   = npy_transport_term_general('B', 64, 16, tstart, tf, B_path_enst)
        np.save(data_path + "B-enst.npy", B)
        del B
        #----------------------------------------------------------------------#
        # Extracting D spectral                                                #
        #----------------------------------------------------------------------#
        D   = npy_transport_term_general('D', 64, 16, tstart, tf, D_path_enst)
        np.save(data_path + "D-enst.npy", D)
        del D
        #----------------------------------------------------------------------#
        # Extracting Pi spectral                                               #
        #----------------------------------------------------------------------#
        Pi   = npy_transport_term_general('Pi', 64, 16, tstart, tf, Pi_path_enst)
        np.save(data_path + "Pi-enst.npy", Pi)
        del Pi
        #----------------------------------------------------------------------#
        # Extracting P spectral                                                #
        #----------------------------------------------------------------------#
        P   = npy_transport_term_general('P', 64, 16, tstart, tf, P_path_enst)
        np.save(data_path + "P-enst.npy", P)
        del P
    #==========================================================================#
    # KE transport terms                                                       #
    #==========================================================================#
    if ke_flag is True:
        #----------------------------------------------------------------------#
        # Extracting A spectral                                                #
        #----------------------------------------------------------------------#
        A   = npy_transport_term_general('A', 64, 16, tstart, tf, A_path_ke)
        np.save(data_path + "A-ke.npy", A)
        del A
        #----------------------------------------------------------------------#
        # Extracting B spectral                                                #
        #----------------------------------------------------------------------#
        B   = npy_transport_term_general('B', 64, 16, tstart, tf, B_path_ke)
        np.save(data_path + "B-ke.npy", B)
        del B
        #----------------------------------------------------------------------#
        # Extracting D spectral                                                #
        #----------------------------------------------------------------------#
        D   = npy_transport_term_general('D', 64, 16, tstart, tf, D_path_ke)
        np.save(data_path + "D-ke.npy", D)
        del D
        #----------------------------------------------------------------------#
        # Extracting Pi spectral                                               #
        #----------------------------------------------------------------------#
        C   = npy_transport_term_general('C', 64, 16, tstart, tf, C_path_ke)
        np.save(data_path + "C-ke.npy", C)
        del C
        #----------------------------------------------------------------------#
        # Extracting P spectral                                                #
        #----------------------------------------------------------------------#
        P   = npy_transport_term_general('P', 64, 16, tstart, tf, P_path_ke)
        np.save(data_path + "P-ke.npy", P)
        del P
    #==========================================================================#
    # Enstrophy and kinetic energy fields                                      #
    #==========================================================================#
    if fields_flag is True:
        #----------------------------------------------------------------------#
        # KE fields                                                            #
        #----------------------------------------------------------------------#
        ke      = npy_transport_term_general('ke', 64, 16, tstart, tf, ke_path)
        np.save(data_path + "ke.npy", ke)
        #----------------------------------------------------------------------#
        # Average kinetic energy                                               #
        #----------------------------------------------------------------------#
        if avg_flag is True:
            ke_avg  = np.loadtxt("../ke-avg.txt", skiprows=2)
            ke_avg  = ke_avg[:,2]
            np.save(data_path + "ke-avg.npy", ke_avg)
            del ke_avg
        del ke
        #----------------------------------------------------------------------#
        # Enstrophy fields                                                     #
        #----------------------------------------------------------------------#
        enst    = npy_transport_term_general('enstrophy', 64, 16, tstart, tf, enst_path)
        np.save(data_path + "enst.npy", enst)
        del enst
        #----------------------------------------------------------------------#
        # Average kinetic energy                                               #
        #----------------------------------------------------------------------#
        if avg_flag is True:
            enst        = npy_transport_term_general('enstrophy', 64, 16, 0, tf, enst_path)
            enst_avg    = ke_average2(enst)
            np.save(data_path + "enst-avg.npy", enst_avg)
            del enst_avg
            del enst

    print("**** Successful Run ****")
    sys.exit(0)

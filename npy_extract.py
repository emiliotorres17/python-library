#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of these subroutines is to extract and compile the data
    from the ALES simulations.

    **** Note: This subroutine has not been unit tested, however a visual
    test has been ran and tested in:
        ~/Project/Fall-2019/blowup/test-10-2019/data-65-8

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
import scipy.io as sci
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Time string                                                             #
#-------------------------------------------------------------------------#
def time_string(
        time):

    """ Generating the time file extension for each file name """
    #---------------------------------------------------------------------#
    # Creating time string                                                #
    #---------------------------------------------------------------------#
    if time < 10:
        time_str    = '_00' + str(time)
    elif time >= 10 and time < 100:
        time_str    = '_0' + str(time)
    else:
        time_str    = '_' + str(time)

    return time_str
#-------------------------------------------------------------------------#
# Processor string                                                        #
#-------------------------------------------------------------------------#
def proc_string(
        proc):                      # processor number

    """ Generating the processor file extension for each file name """
    #---------------------------------------------------------------------#
    # Creating processor string                                           #
    #---------------------------------------------------------------------#
    if proc < 10:
        proc_str    = '_00' + str(proc)
    elif proc >= 10 and proc < 100:
        proc_str    = '_0' + str(proc)

    return proc_str
#-------------------------------------------------------------------------#
# NPY time extraction                                                     #
#-------------------------------------------------------------------------#
def npy_time(
        Nt,                         # number of time steps
        location):                  # file locations

    """ Extracting the time vector from the ALES simulations """
    #---------------------------------------------------------------------#
    # Checking the location path to make sure it has a seperator          #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    #---------------------------------------------------------------------#
    # Extracting the time vector                                          #
    #---------------------------------------------------------------------#
    time        = np.zeros(Nt + 1)
    print_count = 0
    for i in range(0, Nt+1):
        time_str        = time_string(i)
        file_name       = location + 'SimulationTime' + time_str + '.npy'
        time[i]         = np.load(file_name)
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('time extract ---> t_step = %i'       %(i))
            print_count = 0
        print_count += 1

    return time
#-------------------------------------------------------------------------#
# NPY velocity extraction w/ proc = 1                                     #
#-------------------------------------------------------------------------#
def npy_velocity_proc_1(
        ucomp,                          # velocity component
        Nt,                             # number of time steps
        location):                      # data location

    """ Extracting velocity from a simulation done on 1 processor """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep     = os.sep
    proc    = '_000'
    if location[-1] != sep:
        location    = location + sep
    if isinstance(ucomp, str) is False:
        ucomp   = str(ucomp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    u1          = np.zeros((64, 64, 64, Nt+1))
    print_count = 0
    #---------------------------------------------------------------------#
    # Extracting the u-1 data                                             #
    #---------------------------------------------------------------------#
    for i in range(0, Nt+1):
        time_str        = time_string(i)
        file_name       = location + 'Velocity' + ucomp + time_str +\
                            proc + '.npy'
        u1[:,:,:,i]     = np.load(file_name)
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('Velocity' + ucomp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return u1
#-------------------------------------------------------------------------#
# NPY velocity extraction w/ proc = 16                                    #
#-------------------------------------------------------------------------#
def npy_velocity_proc_16(
        ucomp,                          # velocity component
        Nt,                             # number of time steps
        location):                      # data location

    """ Extracting the velocities from a simulation with 16 processors """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(ucomp, str) is False:
        ucomp   = str(ucomp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    u       = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for i in range(0, Nt+1):
        time_str        = time_string(i)
        for n in range(0,16):
            proc        = proc_string(n)
            utemp       = np.load(location + 'Velocity' + ucomp + time_str +\
                                proc + '.npy')
            index       = int(4*n)
            u[index:index + 4, :, :, i]   = utemp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('Velocity' + ucomp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return u
#-------------------------------------------------------------------------#
# Subroutine for extracting the tau (subgrid stress)                      #
#-------------------------------------------------------------------------#
def npy_tau_extract_16(
        tau_comp,                   # subgrid stress component
        Nt,                         # number of time steps
        location):                  # tau file location

    """ Subroutine to extract the subgrib stress """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(tau_comp, str) is False:
        tau_comp   = str(tau_comp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    tau     = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for i in range(0, Nt+1):
        time_str        = time_string(i)
        for n in range(0,16):
            proc        = proc_string(n)
            tau_temp    = np.load(location + 'tau' + tau_comp + time_str +\
                                proc + '.npy')
            index       = int(4*n)
            tau[index:index + 4, :, :, i]   = tau_temp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('tau' + tau_comp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return tau
#-------------------------------------------------------------------------#
# Subroutine for extracting the pressure                                  #
#-------------------------------------------------------------------------#
def npy_pressure_extract_16(
        pressure_comp,              # pressure component
        Nt,                         # number of time steps
        location):                  # location of files

    """ Subroutine to extract the pressure stress """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(pressure_comp, str) is False:
        pressure_comp   = str(pressure_comp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    pressure    = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for i in range(0, Nt+1):
        time_str        = time_string(i)
        for n in range(0,16):
            proc            = proc_string(n)
            pressure_temp   = np.load(location + 'Pressure' + pressure_comp + time_str +\
                                proc + '.npy')
            index           = int(4*n)
            pressure[index:index + 4, :, :, i]   = pressure_temp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('pressure' + pressure_comp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return pressure
#-------------------------------------------------------------------------#
# Subroutine for full extraction and storing time and velocities (npy/mat)#
#-------------------------------------------------------------------------#
def npy_total_extract(
        Nt,                         # number of time steps
        Np,                         # number of processors
        time_path,                  # time file locations
        vel1_path,                  # velocity-1 file location
        vel2_path,                  # velocity-2 file location
        vel3_path,                  # velocity-3 file location
        store_path,                 # storing location
        save_mat    = False):       # flag to save .mat files


    """ Subroutine for extracting and storing time and velocities """
    #---------------------------------------------------------------------#
    # Extracting time                                                     #
    #---------------------------------------------------------------------#
    time    = npy_time(Nt, time_path)
    #---------------------------------------------------------------------#
    # Extracting velocities                                               #
    #---------------------------------------------------------------------#
    if Np == 1:
        u1  = npy_velocity_proc_1(1, Nt, vel1_path)
        u2  = npy_velocity_proc_1(2, Nt, vel2_path)
        u3  = npy_velocity_proc_1(3, Nt, vel3_path)
    elif Np == 16:
        u1  = npy_velocity_proc_16(1, Nt, vel1_path)
        u2  = npy_velocity_proc_16(2, Nt, vel2_path)
        u3  = npy_velocity_proc_16(3, Nt, vel3_path)
    else:
        print("**** Number of processor error \t np = %i ****"      %(np))
    #---------------------------------------------------------------------#
    # Storing .npy files                                                  #
    #---------------------------------------------------------------------#
    np.save(store_path + 'time.npy', time)
    np.save(store_path + 'u1.npy', u1)
    np.save(store_path + 'u2.npy', u2)
    np.save(store_path + 'u3.npy', u3)
    print("\n\n**** Finished storing .npy ****")
    #---------------------------------------------------------------------#
    # Storing .mat files                                                  #
    #---------------------------------------------------------------------#
    if save_mat is not False:
        sci.savemat(store_path + 'time.mat', {'time':time})
        sci.savemat(store_path + 'u1.mat', {'u1':u1})
        sci.savemat(store_path + 'u2.mat', {'u2':u2})
        sci.savemat(store_path + 'u3.mat', {'u3':u3})
        print("\n\n**** Finished storing .mat ****")

    return time, u1, u2, u3
#-------------------------------------------------------------------------#
# NPY velocity extraction w/ proc = 16 from t=t_{0}-t_{f}                 #
#-------------------------------------------------------------------------#
def npy_velocity_interval(
        ucomp,                          # velocity component
        tstart,                         # time start (i.e., 100)
        tfinish,                        # time finish (i.e., 500)
        location):                      # data location

    """ Extracting the velocities from a simulation with 16 processors
    from time range t=t_{0} to t=t_{finish}"""
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(ucomp, str) is False:
        ucomp   = str(ucomp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    Nt      = int(tfinish-tstart)
    u       = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for count, i in enumerate(range(tstart, tfinish+1)):
        time_str        = time_string(i)
        for n in range(0,16):
            proc        = proc_string(n)
            utemp       = np.load(location + 'Velocity' + ucomp + time_str +\
                                proc + '.npy')
            index       = int(4*n)
            u[index:index + 4, :, :, count]   = utemp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('Velocity' + ucomp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return u
#-------------------------------------------------------------------------#
# NPY time extraction for a defined time interval                         #
#-------------------------------------------------------------------------#
def npy_time_interval(
        t0,                         # start time
        tf,                         # finish time
        location):                  # file locations

    """ Extracting the time vector from the ALES simulations for a given
    time interval """
    #---------------------------------------------------------------------#
    # Checking the location path to make sure it has a seperator          #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    #---------------------------------------------------------------------#
    # Extracting the time vector                                          #
    #---------------------------------------------------------------------#
    Nt          = int(tf-t0)
    time        = np.zeros(Nt + 1)
    print_count = 0
    for count, i in enumerate(range(t0, tf+1)):
        time_str        = time_string(i)
        file_name       = location + 'SimulationTime' + time_str + '.npy'
        time[count]     = np.load(file_name)
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('time extract ---> t_step = %i'       %(i))
            print_count = 0
        print_count += 1

    return time
#-------------------------------------------------------------------------#
# Subroutine for extracting the pressure for a certain time interval      #
#-------------------------------------------------------------------------#
def npy_pressure_interval(
        pressure_comp,              # pressure component
        tstart,                     # time start
        tfinish,                    # finish time
        location):                  # location of files

    """ Subroutine to extract the pressure for a given time interval """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(pressure_comp, str) is False:
        pressure_comp   = str(pressure_comp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    Nt          = int(tfinish-tstart)
    pressure    = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for count, i in enumerate(range(tstart, tfinish+1)):
        time_str        = time_string(i)
        for n in range(0,16):
            proc            = proc_string(n)
            pressure_temp   = np.load(location + 'Pressure' + pressure_comp + time_str +\
                                proc + '.npy')
            index           = int(4*n)
            pressure[index:index + 4, :, :, count]   = pressure_temp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('pressure' + pressure_comp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return pressure
#-------------------------------------------------------------------------#
# Subroutine for extracting the tau for a given time interval             #
#-------------------------------------------------------------------------#
def npy_tau_interval(
        tau_comp,                   # subgrid stress component
        tstart,                     # starting time
        tfinish,                    # finishing time
        location):                  # tau file location

    """ Subroutine to extract the subgrib stress for a given time
    interval """
    #---------------------------------------------------------------------#
    # Defining file variables                                             #
    #---------------------------------------------------------------------#
    sep = os.sep
    if location[-1] != sep:
        location    = location + sep
    if isinstance(tau_comp, str) is False:
        tau_comp   = str(tau_comp)
    #---------------------------------------------------------------------#
    # Preallocating space for the velocity                                #
    #---------------------------------------------------------------------#
    Nt      = int(tfinish-tstart)
    tau     = np.zeros((64, 64, 64, Nt+1))
    #---------------------------------------------------------------------#
    # Extracting the data                                                 #
    #---------------------------------------------------------------------#
    print_count = 0
    for count, i in enumerate(range(tstart, tfinish+1)):
        time_str        = time_string(i)
        for n in range(0,16):
            proc        = proc_string(n)
            tau_temp    = np.load(location + 'tau' + tau_comp + time_str +\
                                proc + '.npy')
            index       = int(4*n)
            tau[index:index + 4, :, :, count]   = tau_temp
        #-----------------------------------------------------------------#
        # Printing time step output                                       #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('tau' + tau_comp + ' ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return tau
#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == "__main__":

    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(["clear"])
    #---------------------------------------------------------------------#
    # Print statement                                                     #
    #---------------------------------------------------------------------#
    print('Has not been united tested')

    sys.exit(0)

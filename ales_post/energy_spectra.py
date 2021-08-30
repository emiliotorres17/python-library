#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to find the energy spectra of the kinetic 
    energy for ALES simulation. 
    
    **** Note: Currently hard coded for N = 64 ****
    
    **** Has not been tested yet! ****

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Purpose                                                                 #
#=========================================================================#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import os
import sys
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from ales_post.plot_settings        import plot_setting
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Shell average                                                           #
#-------------------------------------------------------------------------#
def shell_average(E3, km):
    """
    Compute the 1D, shell-averaged, spectrum of the 3D Fourier-space
    variable E3.

    Arguments:
        comm - MPI intracommunicator
        nk   - scalar length of 1-D spectrum
        km   - wavemode of each n-D wavevector
        E3   - 3-dimensional complex or real Fourier-space scalar
    """
    nz, nny, nk = E3.shape
    E1 = zeros(nk, dtype=E3.dtype)
    z = zeros_like(E3)

    if km[0, 0, 0] == 0:
        # this mpi task has the DC mode, only works for 1D domain decomp
        E1[0] = E3[0, 0, 0]

    for k in range(1, nk):
        E1[k] = psum(where(km==k, E3, z))

    #comm.Allreduce(MPI.IN_PLACE, E1, op=MPI.SUM)

    return E1
#-------------------------------------------------------------------------#
# 3D real-valued FFTs                                                     #
#-------------------------------------------------------------------------#
def rfft3(u, fu=None):
    """
    Compute MPI-distributed, real-to-complex 3D FFT.
    Input array must have only three dimensions (not curently checked)
    temp and fu are complex data work arrays where the array view for fu can
    be passed in from the calling function
    """
    ntasks = 1
    nnz, ny, nx = u.shape
    nk = nx//2+1
    nny = ny//ntasks
    nz = nnz*ntasks

    if fu is None:
        fu = empty([nz, nny, nk], dtype=complex128)

    temp1 = empty([nnz, ny, nk], dtype=complex128)

    temp1[:] = fft.rfft2(u, axes=(1, 2))
    fu[:] = rollaxis(temp1.reshape([nnz, ntasks, nny, nk]),
                        1).reshape(fu.shape)
    fu[:] = fft.fft(fu, axis=0)

    return fu
#-------------------------------------------------------------------------#
# Energy spectra                                                          #
#-------------------------------------------------------------------------#
def energy_spectra(
        u1,
        u2,
        u3):
    
    """ Calculating the energy spectra """
    #---------------------------------------------------------------------#
    # Defining size of filter variables                                   #
    #---------------------------------------------------------------------#
    nnk     = zeros(3, dtype=int)
    nnk[0]  = 64
    nnk[1]  = 64
    nnk[2]  = 33
    U_hat   = empty((3, *nnk), dtype=complex)  # solution vector
    #---------------------------------------------------------------------#
    # Taking FFT                                                          #
    #---------------------------------------------------------------------#
    rfft3(u1, U_hat[0])
    rfft3(u2, U_hat[1])
    rfft3(u3, U_hat[2])
    #---------------------------------------------------------------------#
    # Defining variables to define the K vector                           #
    #---------------------------------------------------------------------#
    iks     = zeros(3, dtype=int)
    ike     = iks+nnk
    nx      = zeros(3, dtype=int)
    nx[0]   = 64
    nx[1]   = 64
    nx[2]   = 64
    #---------------------------------------------------------------------#
    # Defining K vector                                                   #
    #---------------------------------------------------------------------#
    k0      = fft.fftfreq(nx[0])*nx[0]
    k1      = fft.fftfreq(nx[1])[iks[1]:ike[1]]*nx[1]
    k2      = fft.rfftfreq(nx[2])*nx[2]
    K       = array(meshgrid(k0, k1, k2, indexing='ij'), dtype=int)
    Ksq     = sum(square(K), axis=0)
    Kmod    = floor(sqrt(Ksq)).astype(int)
    with errstate(divide='ignore'):
        K_Ksq = where(Ksq > 0, K/Ksq, 0.0)
    #---------------------------------------------------------------------#
    # Output kinetic energy spectrum                                      #
    #---------------------------------------------------------------------#
    spect3d             = sum(real(U_hat*conj(U_hat)), axis=0)
    spect3d[..., 0]     *= 0.5
    spect1d             = shell_average(spect3d, Kmod)
    #---------------------------------------------------------------------#
    # Defining corresponding wave numbers                                 #
    #---------------------------------------------------------------------#
    k   = linspace(1, 33, 33)

    return k, spect1d
#-------------------------------------------------------------------------#
# Energy spectra plot                                                     #
#-------------------------------------------------------------------------#
def plot_energy_spectra(
        u1,                   
        u2,                   
        u3,
        symbol,
        Lab,
        settings_flag = False,
        name = False):

    """ Plotting the energy spectra """
    #---------------------------------------------------------------------#
    # Name settings                                                       #
    #---------------------------------------------------------------------#
    if settings_flag is False and name is False:
        name = 'energy-spectra.png'
    #---------------------------------------------------------------------#
    # Finding the energy spectrum                                         #
    #---------------------------------------------------------------------#
    [k, sepct1d]    = energy_spectra(u1, u2, u3)
    #---------------------------------------------------------------------#
    # Plotting spectra                                                    #
    #---------------------------------------------------------------------#
    plot_setting()
    plt.loglog(k[1:32], spect1d[1:32]/spect1d[1], symbol, lw=1.5, label= Lab) 
    #---------------------------------------------------------------------#
    # Plot settings                                                       #
    #---------------------------------------------------------------------#
    if settings_flag is not False:
        plt.ylim([10**-3, 2])
        plt.xlim([1, 60])
        plt.ylabel('$E(k)/E_{0}$')
        plt.xlabel('$k$')
        plt.grid(True, which='both', ls='--', c='gray')
        plt.legend(loc=0)
        plt.savefig(name)

    return

#=========================================================================#
# Main                                                                    #
#=========================================================================#
if __name__ == '__main__':
    #---------------------------------------------------------------------#
    # Main preamble                                                       #
    #---------------------------------------------------------------------#
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = '/home/eetorres/MAE-792/Spring-2021/june/tf5/TF5_for_Emilio/TF5-CS-0.85-no-limiting/data/'
    media_path  = pwd + '%c..%cmedia%c'             %(sep, sep, sep)
    #---------------------------------------------------------------------#
    # load data                                                           #
    #---------------------------------------------------------------------#
    u1          = load(data_path + 'velocity1.npy')[:,:,:,-1]
    u2          = load(data_path + 'velocity2.npy')[:,:,:,-1]
    u3          = load(data_path + 'velocity3.npy')[:,:,:,-1]

    bs  = scipy.io.loadmat('BS-spectra.mat')
    bs  = bs['Ek_85']
    kbs = scipy.io.loadmat('k-bs.mat')
    kbs = transpose(kbs['k'])




    print('**** Successful run ****')
    sys.exit(0)

#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this script is to calculate the non-local and local
    strain rates from the total strain rate field using the methods in Self
    Attenuation (spectral) and PRE (gradient) papers.

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
# Spectral methods                                                        #
#-------------------------------------------------------------------------#
def spectral_non_local(
        s,
        rad):
    """ 
        Finding the non-local and local strain rates from the total strain
        using the method in the Self Attenuation paper's method
    """
    #---------------------------------------------------------------------#
    # Generating the wave number                                          #
    #---------------------------------------------------------------------#
    dim     = s.shape[0] 
    k       = np.fft.fftfreq(dim) * dim
    Kfield  = np.array(np.meshgrid(k,k,k,indexing='ij'))
    kmag    = np.sqrt(np.sum(np.square(Kfield.astype(np.float64)), axis=0))
    kmag[0,0,0] = 1.0
    #---------------------------------------------------------------------#
    # Calculating the Fourier transform terms                             #
    #---------------------------------------------------------------------#
    f       = (3.0*(np.sin(kmag*rad)-kmag*rad*np.cos(kmag*rad)))\
                    /((kmag*rad)**3.0)
    sNL     = np.fft.ifftn(f*np.fft.fftn(s)).real 
    sL      = s - sNL

    return (sNL, sL)
#-------------------------------------------------------------------------#
# Gradient methods                                                        #
#-------------------------------------------------------------------------#
def gradient_non_local(
        s,
        rad):
    """ 
        Finding the non-local and local strain rates from the total strain
        using the method in the PRE paper 
    """
    #---------------------------------------------------------------------#
    # Setting up the wave number for differentiation                      #
    #---------------------------------------------------------------------#
    dim         = s.shape[0]
    kspec       = np.fft.fftfreq(dim) * dim
    Kfield      = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
    Ksq         = np.sum(np.square(Kfield), axis=0)
    #---------------------------------------------------------------------#
    # Finding the second term of the non-local strain                     #
    #---------------------------------------------------------------------#
    s_hat       = np.fft.fftn(s)
    term2_hat   = -Ksq*s_hat
    term2       = np.fft.ifftn(term2_hat).real
    #---------------------------------------------------------------------#
    # Finding the third term of the non-local strain                      #
    #---------------------------------------------------------------------#
    term3_1     = np.fft.ifftn(-Ksq*np.fft.fftn(s)).real
    term3_hat   = -Ksq*np.fft.fftn(term3_1)
    term3       = np.fft.ifftn(term3_hat).real
    #---------------------------------------------------------------------#
    # Multiplying both terms by their respective constants                #
    #---------------------------------------------------------------------#
    term2   *= R**2.0/10.0
    term3   *= R**4.0/280.0
    #---------------------------------------------------------------------#
    # Finding the no-local and local strain                               #
    #---------------------------------------------------------------------#
    sNL     = s + term2 + term3
    sL      = s - sNL

    return (sNL, sL)
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
    data_path   = pwd + '%c..%cdata%c'          %(sep, sep, sep)
    media_path  = pwd + '%c..%cmedia%c'         %(sep, sep, sep)
    
    print('**** This has not been unit tested ****')

    print('**** Successful run ****')
    sys.exit(0)

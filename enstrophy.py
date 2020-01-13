#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to calculate the enstrophy field.

Author:
    Emilio Torres
========================================================================"""
#=========================================================================#
# Preamble                                                                #
#=========================================================================#
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import sys
from subprocess import call
import numpy as np
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Enstrophy                                                               #
#-------------------------------------------------------------------------#
def enstrophy_field(
        omega1,                     # vorticity-1 component
        omega2,                     # vorticity-2 component
        omega3):                    # vorticity-3 component

    """ Calculating the enstrophy field """
    #---------------------------------------------------------------------#
    # Determining the shape of the enstrophy field                        #
    #---------------------------------------------------------------------#
    dim     = omega1.shape
    ens     = np.empty((dim[0], dim[1], dim[2], dim[3]))
    #---------------------------------------------------------------------#
    # Time loop                                                           #
    #---------------------------------------------------------------------#
    print_count = 0
    for i in range(0, dim[3]):
        term1           = np.square(omega1[:,:,:,i])
        term2           = np.square(omega2[:,:,:,i])
        term3           = np.square(omega3[:,:,:,i])
        ens[:,:,:,i]    = 0.5*(term1 + term2 + term3)
        #-----------------------------------------------------------------#
        # Printing statement                                              #
        #-----------------------------------------------------------------#
        if print_count > 20:
            print('Enstrophy field ---> t_step = %i'      %(i))
            print_count = 0
        print_count += 1

    return ens
#-------------------------------------------------------------------------#
# Q term from Dahm's emailed notes                                        #
#-------------------------------------------------------------------------#
def Q_term(
        omega1,             # vorticity-1 component
        omega2,             # vorticity-2 component
        omega3,             # vorticity-3 component
        s11,                # strain rate-11 component
        s12,                # strain rate-12 component
        s13,                # strain rate-13 component
        s22,                # strain rate-22 component
        s23,                # strain rate-23 component
        s33):               # strain rate-33 component

    """ Calculating the Q term from Dahm's note """
    #---------------------------------------------------------------------#
    # Numerator and denominator                                           #
    #---------------------------------------------------------------------#
    num     = omega1*s11*omega1 + omega1*s12*omega2 +  omega1*s13*omega3 +\
                omega2*s12*omega1 + omega2*s22*omega2 + omega2*s23*omega3+\
                omega3*s13*omega1 + omega3*s23*omega2 + omega3*s33*omega3
    den1    = omega1*omega1 + omega2*omega2 + omega3*omega3
    den2    = (s11*s11 + s12*s12 + s13*s13 + s12*s12 + s22*s22 + s23*s23 +\
                s13*s13 + s23*s23 + s33*s33)**0.5
    den     = ((2.0/3.0)**0.5)* den1 * den2
    #---------------------------------------------------------------------#
    # Q calculation                                                       #
    #---------------------------------------------------------------------#
    Q       = num/den

    return Q
#-------------------------------------------------------------------------#
# R term from Dahm's emailed notes                                        #
#-------------------------------------------------------------------------#
def R_term(
        enst,               # enstrophy field
        omega1,             # vorticity-1 component
        omega2,             # vorticity-2 component
        omega3,             # vorticity-3 component
        s11,                # strain rate-11 component
        s12,                # strain rate-12 component
        s13,                # strain rate-13 component
        s22,                # strain rate-22 component
        s23,                # strain rate-23 component
        s33,                # strain rate-33 component
        diff = False):      # differentiation flag

    """ Calculating the R term from Dahm's scanned notes
            **** Note: need to update the gradient terms to use spectral
            methods """
    #---------------------------------------------------------------------#
    # Defining domain variables                                           #
    #---------------------------------------------------------------------#
    pi      = np.pi                             # pi
    dx      = (2.0*pi)/64.0                     # spatial step
    nu      = 0.000185                          # default viscosity
    #---------------------------------------------------------------------#
    # Spectral differentiation variables                                  #
    #---------------------------------------------------------------------#
    dim     = 64
    kspec   = np.fft.fftfreq(dim) * dim
    Kfield  = np.array(np.meshgrid(kspec, kspec, kspec, indexing='ij'))
    #---------------------------------------------------------------------#
    # Spectral differentiation variables                                  #
    #---------------------------------------------------------------------#
    term1   = np.zeros((dim, dim, dim))
    term2   = np.zeros((dim, dim, dim))
    term3   = np.zeros((dim, dim, dim))
    #---------------------------------------------------------------------#
    # Numerator (numpy gradient tool)                                     #
    #---------------------------------------------------------------------#
    if diff is not False:
        term1   = np.gradient(enst,dx, edge_order=2)[0]
        term2   = np.gradient(enst,dx, edge_order=2)[1]
        term3   = np.gradient(enst,dx, edge_order=2)[2]
    #---------------------------------------------------------------------#
    # Numerator (spectral differentiation)                                #
    #---------------------------------------------------------------------#
    else:
        term1   = 0.5*np.fft.ifftn(1j*Kfield[0]*np.fft.fftn(enst) +\
                        1j*Kfield[0]*np.fft.fftn(enst)).real
        term2   = 0.5*np.fft.ifftn(1j*Kfield[1]*np.fft.fftn(enst) +\
                        1j*Kfield[1]*np.fft.fftn(enst)).real
        term3   = 0.5*np.fft.ifftn(1j*Kfield[2]*np.fft.fftn(enst) +\
                        1j*Kfield[2]*np.fft.fftn(enst)).real
    #---------------------------------------------------------------------#
    # Numerator                                                           #
    #---------------------------------------------------------------------#
    num     = nu*(term1**2.0+ term2**2.0 + term3**2.0)
    #---------------------------------------------------------------------#
    # Denominator                                                         #
    #---------------------------------------------------------------------#
    den     = omega1*s11*omega1 + omega1*s12*omega2 +  omega1*s13*omega3 +\
                omega2*s12*omega1 + omega2*s22*omega2 + omega2*s23*omega3+\
                omega3*s13*omega1 + omega3*s23*omega2 + omega3*s33*omega3
    #---------------------------------------------------------------------#
    # R calculation                                                       #
    #---------------------------------------------------------------------#
    R       = num/den

    return R
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
    print('This has not been unit tested')

    sys.exit(0)

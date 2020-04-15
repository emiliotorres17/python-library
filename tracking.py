#!/usr/bin/env python3
"""=============================================================================
Purpose:
    The purpose of this subroutine is to perform the tracking of a particle
    from an ALES simulation.

    **** Warning: This has not been tested, and there is no proof it works as
                  intended.
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
from  min_max       import find_max3D
from  min_max       import find_min3D
#from normalize      import norm_min_max
#==============================================================================#
# User defined functions                                                       #
#==============================================================================#
#------------------------------------------------------------------------------#
# Generating color fields                                                      #
#------------------------------------------------------------------------------#
def contour_gen(
        field,
        Xin,
        Yin,
        Zin,
        tval):

    """ Generating the contours in order to visually track """
    #--------------------------------------------------------------------------#
    # domain variables                                                         #
    #--------------------------------------------------------------------------#
    print(Xin)
    print(Yin)
    print(Zin)
    xvec        = np.linspace(0.0, 2.0*np.pi, 64)
    [X1, X2]    = np.meshgrid(xvec,xvec)
    field       = field[:,:,:,tval]
    #--------------------------------------------------------------------------#
    # colorbar settings                                                        #
    #--------------------------------------------------------------------------#
    Vmin        = np.amin(field[:,:,Zin])
    Vmax        = np.amax(field[:,:,Zin])
    print(field[Xin, Yin, Zin])
    print(Vmax)
    #sys.exit(100)
    Dp          = (Vmax - Vmin)/500.0
    #--------------------------------------------------------------------------#
    # generating plots                                                         #
    #--------------------------------------------------------------------------#
    fig0, ax0   = plt.subplots(1,1,)
    cnt         = ax0.contourf(X1, X2, field[:,:,Zin], np.arange(Vmin, Vmax, Dp),\
                            cmap='jet', extend='both')
    for c in cnt.collections:
        c.set_edgecolors("face")
    cbar0   = plt.colorbar(cnt)
    #--------------------------------------------------------------------------#
    # saving picture                                                           #
    #--------------------------------------------------------------------------#
    plt.savefig('media/Dk_DT-' + str(int(tval)) + '.png')
    plt.clf()

    return
#------------------------------------------------------------------------------#
# Periodic boundary conditions                                                 #
#------------------------------------------------------------------------------#
def location_adjustment(
        loc):

    """ Subroutine to account for the periodic boundary conditions """
    #--------------------------------------------------------------------------#
    # Domain variables                                                         #
    #--------------------------------------------------------------------------#
    PI  = np.pi
    #--------------------------------------------------------------------------#
    # Location >  2pi                                                          #
    #--------------------------------------------------------------------------#
    if loc > 2.0*PI:
        Loc = loc - 2.0*PI
    #--------------------------------------------------------------------------#
    # Location <  0                                                            #
    #--------------------------------------------------------------------------#
    elif loc < 0.0:
        Loc = loc + 2.0*PI
    #--------------------------------------------------------------------------#
    # 0.0 <= Location <= 2.0*pi                                                #
    #--------------------------------------------------------------------------#
    else:
        Loc = loc

    return Loc
#------------------------------------------------------------------------------#
# Subroutine to find specific value                                            #
#------------------------------------------------------------------------------#
def find_value(
        xvec,                   # x vector
        xval):                  # x value

    """ subroutine to find the coordinate value """
    #--------------------------------------------------------------------------#
    # Finding each distance between each point                                 #
    #--------------------------------------------------------------------------#
    dist     = np.zeros(len(xvec))
    for Count, Val in enumerate(xvec):
        dist[Count]  = abs(Val-xval)
    #--------------------------------------------------------------------------#
    # Finding the minimum distance                                             #
    #--------------------------------------------------------------------------#
    index   = np.unravel_index(np.argmin(dist, axis=None), dist.shape)[0]
    index   = int(index)

    return index
#------------------------------------------------------------------------------#
# Tracking subroutine                                                          #
#------------------------------------------------------------------------------#
def tracking(
        start,                      # start time point
        end,                        # end time point
        Time,                       # time vector
        loc_vec,                    # x, y, and z vectors
        crd_vec,                    # x, y, and z starting coordinates
        Ux,                         # x-velocity
        Uy,                         # y-velocity
        Uz,                         # z-velocity
        BA_flag = True,             # box average flag
        debug   = False):           # debug flag

    """ Subroutine to track a particular point """
    #--------------------------------------------------------------------------#
    # Setting the x, y, and z vectors and coordinates                          #
    #--------------------------------------------------------------------------#
    Xvec    = loc_vec[0]
    Yvec    = loc_vec[1]
    Zvec    = loc_vec[2]
    Nx      = crd_vec[0]
    Ny      = crd_vec[1]
    Nz      = crd_vec[2]
    #--------------------------------------------------------------------------#
    # Initial conditions                                                       #
    #--------------------------------------------------------------------------#
    xt  = Xvec[Nx]
    yt  = Yvec[Ny]
    zt  = Zvec[Nz]
    #--------------------------------------------------------------------------#
    # Preallocating time loop variables                                        #
    #--------------------------------------------------------------------------#
    x_pts   = np.zeros((end-start)+1)
    y_pts   = np.zeros((end-start)+1)
    z_pts   = np.zeros((end-start)+1)
    #--------------------------------------------------------------------------#
    # Time loop                                                                #
    #--------------------------------------------------------------------------#
    for t in range(end, start-1, -1):
        #----------------------------------------------------------------------#
        # Storing coordinate values                                            #
        #----------------------------------------------------------------------#
        x_pts[t-1]  = Nx
        y_pts[t-1]  = Ny
        z_pts[t-1]  = Nz
        #----------------------------------------------------------------------#
        # previous x, y, and z values                                          #
        #----------------------------------------------------------------------#
        xt_old  = xt
        yt_old  = yt
        zt_old  = zt
        #----------------------------------------------------------------------#
        # Previous Nx, Ny, Nz                                                  #
        #----------------------------------------------------------------------#
        Nx_old  = Nx
        Ny_old  = Ny
        Nz_old  = Nz
        #----------------------------------------------------------------------#
        # Finding x, y, and z velocities                                       #
        #----------------------------------------------------------------------#
        if BA_flag is True:
            vel_x   = box_average(Ux[:,:,:,t],Nx, Ny, Nz)   # corresponding x-vel.
            vel_y   = box_average(Uy[:,:,:,t],Nx, Ny, Nz)   # corresponding x-vel.
            vel_z   = box_average(Uz[:,:,:,t],Nx, Ny, Nz)   # corresponding x-vel.
        else:
            vel_x   = star_average(Ux[:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
            vel_y   = star_average(Uy[:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
            vel_z   = star_average(Uz[:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
        #----------------------------------------------------------------------#
        # Calculating new locations                                            #
        #----------------------------------------------------------------------#
        dt      = Time[t] - Time[t-1]       # Delta t
        xt      = xt + vel_x*dt             # x location at t^{n-1}
        yt      = yt + vel_y*dt             # y location at t^{n-1}
        zt      = zt + vel_z*dt             # z location at t^{n-1}
        #----------------------------------------------------------------------#
        # Applying periodic boundary conditions                                #
        #----------------------------------------------------------------------#
        xt      = location_adjustment(xt)
        yt      = location_adjustment(yt)
        zt      = location_adjustment(zt)
        #----------------------------------------------------------------------#
        # Finding the x, y, and z coordinate                                   #
        #----------------------------------------------------------------------#
        Nx      = find_value(Xvec, xt)
        Ny      = find_value(Yvec, yt)
        Nz      = find_value(Zvec, zt)
        #----------------------------------------------------------------------#
        # Debug print                                                          #
        #----------------------------------------------------------------------#
        if debug is True:
            print('x-values:')
            print('\ttime=%.5f\tdt=%.5e\tvel_{x}=%.5f\tx_{n+1}=%.5f\tx_{n}=%.5f\tNx_{n+1}=%i\tNx_{n}=%i'\
                    %(t,dt,vel_x,xt_old,xt,Nx_old,Nx))
            print('y-values:')
            print('\ttime=%.5f\tdt=%.5e\tvel_{y}=%.5f\ty_{n+1}=%.5f\ty_{n}=%.5f\tNy_{n+1}=%i\tNy_{n}=%i'\
                    %(t,dt,vel_y,yt_old,yt,Ny_old,Ny))
            print('z-values:')
            print('\ttime=%.5f\tdt=%.5e\tvel_{z}=%.5f\tz_{n+1}=%.5f\tz_{n}=%.5f\tNz_{n+1}=%i\tNz_{n}=%i'\
                    %(t,dt,vel_z,zt_old,zt,Nz_old,Nz))

    return x_pts, y_pts, z_pts
#------------------------------------------------------------------------------#
# Star averaging                                                               #
#------------------------------------------------------------------------------#
def star_average(
        field,
        Nx,
        Ny,
        Nz):

    """ Calculating a star average of a single point in the field """
    #--------------------------------------------------------------------------#
    # Calculating the plus (i+1) points                                        #
    #--------------------------------------------------------------------------#
    Nx_plus     = coordinate_adjustment(Nx+1)
    Ny_plus     = coordinate_adjustment(Ny+1)
    Nz_plus     = coordinate_adjustment(Nz+1)
    #--------------------------------------------------------------------------#
    # Calculating the minus (i-1) points                                       #
    #--------------------------------------------------------------------------#
    Nx_minus   = coordinate_adjustment(Nx-1)
    Ny_minus   = coordinate_adjustment(Ny-1)
    Nz_minus   = coordinate_adjustment(Nz-1)
    #--------------------------------------------------------------------------#
    # Finding the values of seven points                                       #
    #--------------------------------------------------------------------------#
    Vals    = np.zeros(7)
    Vals[0] = field[Nx,Ny,Nz]
    Vals[1] = field[Nx_plus,Ny,Nz]
    Vals[2] = field[Nx_minus,Ny,Nz]
    Vals[3] = field[Nx,Ny_plus,Nz]
    Vals[4] = field[Nx,Ny_minus,Nz]
    Vals[5] = field[Nx,Ny,Nz_plus]
    Vals[6] = field[Nx,Ny,Nz_minus]
    Value   = np.mean(Vals)

    return Value
#------------------------------------------------------------------------------#
# Periodic boundary conditions for coordinates                                 #
#------------------------------------------------------------------------------#
def coordinate_adjustment(
        coor):

    """ Adjusting the coordinate to range from 0 <= coor <= 63 """
    #--------------------------------------------------------------------------#
    # coordinates > 63                                                         #
    #--------------------------------------------------------------------------#
    if coor > 63:
        Coor = coor - 63
    elif coor < 0:
        Coor = 63 + coor
    else:
        Coor = coor

    return Coor
#------------------------------------------------------------------------------#
# Box averaging                                                                #
#------------------------------------------------------------------------------#
def box_average(
        field,                  # field
        Nx,                     # x-coordinate
        Ny,                     # y-coordinate
        Nz):                    # z-coordinate

    """ Calculating the average of a 3x3x3 box around the point of interest """
    #--------------------------------------------------------------------------#
    # i-pts                                                                    #
    #--------------------------------------------------------------------------#
    xpt         = coordinate_adjustment(Nx)
    xpt_p1      = coordinate_adjustment(Nx+1)
    xpt_m1      = coordinate_adjustment(Nx-1)
    #--------------------------------------------------------------------------#
    # j-pts                                                                    #
    #--------------------------------------------------------------------------#
    ypt         = coordinate_adjustment(Ny)
    ypt_p1      = coordinate_adjustment(Ny+1)
    ypt_m1      = coordinate_adjustment(Ny-1)
    #--------------------------------------------------------------------------#
    # k-pts                                                                    #
    #--------------------------------------------------------------------------#
    zpt         = coordinate_adjustment(Nz)
    zpt_p1      = coordinate_adjustment(Nz+1)
    zpt_m1      = coordinate_adjustment(Nz-1)
    #--------------------------------------------------------------------------#
    # Finding the values and mean of the 27 pts                                #
    #--------------------------------------------------------------------------#
    Vals        = np.zeros(27)
    Vals[0]     = field[xpt,    ypt, zpt]
    Vals[1]     = field[xpt_p1, ypt, zpt]
    Vals[2]     = field[xpt_m1, ypt, zpt]
    Vals[3]     = field[xpt,    ypt, zpt_m1]
    Vals[4]     = field[xpt,    ypt, zpt_p1]
    Vals[5]     = field[xpt_m1, ypt, zpt_m1]
    Vals[6]     = field[xpt_p1, ypt, zpt_m1]
    Vals[7]     = field[xpt_m1, ypt, zpt_p1]
    Vals[8]     = field[xpt_p1, ypt, zpt_p1]
    Vals[9]     = field[xpt,    ypt_p1, zpt]
    Vals[10]    = field[xpt_p1, ypt_p1, zpt]
    Vals[11]    = field[xpt_m1, ypt_p1, zpt]
    Vals[12]    = field[xpt,    ypt_p1, zpt_m1]
    Vals[13]    = field[xpt,    ypt_p1, zpt_p1]
    Vals[14]    = field[xpt_m1, ypt_p1, zpt_m1]
    Vals[15]    = field[xpt_p1, ypt_p1, zpt_m1]
    Vals[16]    = field[xpt_m1, ypt_p1, zpt_p1]
    Vals[17]    = field[xpt_p1, ypt_p1, zpt_p1]
    Vals[18]    = field[xpt,    ypt_m1, zpt]
    Vals[19]    = field[xpt_p1, ypt_m1, zpt]
    Vals[20]    = field[xpt_m1, ypt_m1, zpt]
    Vals[21]    = field[xpt,    ypt_m1, zpt_m1]
    Vals[22]    = field[xpt,    ypt_m1, zpt_p1]
    Vals[23]    = field[xpt_m1, ypt_m1, zpt_m1]
    Vals[24]    = field[xpt_p1, ypt_m1, zpt_m1]
    Vals[25]    = field[xpt_m1, ypt_m1, zpt_p1]
    Vals[26]    = field[xpt_p1, ypt_m1, zpt_p1]
    Vals        = np.mean(Vals)

    return Vals
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
    data_path   = 'data%c'                  %(sep)
    #---------------------------------------------------------------------#
    # Loading data                                                        #
    #---------------------------------------------------------------------#
    print('Loading data:')
    time    = np.load(data_path + 'time.npy')
    print('\ttime')
    u_x     = np.load(data_path + 'vel_x.npy')
    print('\tx-veloccity')
    u_y     = np.load(data_path + 'vel_y.npy')
    print('\ty-velocity')
    u_z     = np.load(data_path + 'vel_z.npy')
    print('\tz-velocity')
    Dk      = np.load(data_path + 'Dk_Dt.npy')
    print('\tDk/Dt')
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    x   = np.linspace(0, 2.0*np.pi, 64)
    y   = np.linspace(0, 2.0*np.pi, 64)
    z   = np.linspace(0, 2.0*np.pi, 64)
    #---------------------------------------------------------------------#
    # Finding the maximum value of the field                              #
    #---------------------------------------------------------------------#
    tfinal          = 19
    [val, X, Y, Z]  = find_max3D(Dk[:,:,:,tfinal])
    #---------------------------------------------------------------------#
    # Generating contours plots starting at T=t_{f}                       #
    #---------------------------------------------------------------------#
    for i in range(tfinal, tfinal-6, -1):
        val = find_max3D(Dk[:,:,:,i])[0]
        print('maximum value=%.7f\tx-loc=%.5f\ty-loc=%.5f\tz-loc=%.5f'\
                    %(val, X, Y, Z))
        contour_gen(Dk, X, Y, Z, i)
    #---------------------------------------------------------------------#
    # Tracking                                                            #
    #---------------------------------------------------------------------#
    [X,Y,Z] = tracking(1, tfinal, time, (x,y,z), [X,Y,Z], u_x, u_y,\
                        u_z, True, True)
    print(X)
    print(Y)
    print(Z)

    print('**** Successful run ****')
    sys.exit(0)

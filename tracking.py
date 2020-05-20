#!/usr/bin/env python3
"""========================================================================
Purpose:
    The purpose of this subroutine is to perform the tracking of a particle
    from an ALES simulation.

    **** Warning: This has not been tested, and there is no proof it works as
                  intended.                                         
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
#-------------------------------------------------------------------------#
# User packages                                                           #
#-------------------------------------------------------------------------#
from  min_max       import find_max3D
from  min_max       import find_min3D
#=========================================================================#
# User defined functions                                                  #
#=========================================================================#
#-------------------------------------------------------------------------#
# Distance formula                                                        #
#-------------------------------------------------------------------------#
def distance(
        x1,
        x2,
        y1,
        y2,
        z1,
        z2):

    """ Calculating the distance between two points """
    #---------------------------------------------------------------------#
    # Calculating distance                                                #
    #---------------------------------------------------------------------#
    dist    = np.sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)

    return dist
#-------------------------------------------------------------------------#
# Tracking subroutine 2                                                   #
#-------------------------------------------------------------------------#
def tracking2(
        start,                      # start time point
        end,                        # end time point
        Time,                       # time vector
        loc_vec,                    # x, y, and z vectors
        crd_vec,                    # x, y, and z starting coordinates
        vel,                        # vel_x, vel_y, vel_z
        debug   = False):           # debug flag

    """ Subroutine to track a particle part dos """
    #---------------------------------------------------------------------#
    # Preallocating variables                                             #
    #---------------------------------------------------------------------#
    dist        = np.zeros((64,64,64))
    xpt         = np.zeros(end-start)
    ypt         = np.zeros(end-start)
    zpt         = np.zeros(end-start)
    #---------------------------------------------------------------------#
    # Setting coordinates                                                 #
    #---------------------------------------------------------------------#
    Nx          = crd_vec[0]
    Ny          = crd_vec[1]
    Nz          = crd_vec[2]
    #---------------------------------------------------------------------#
    # setting the target values                                           #
    #---------------------------------------------------------------------#
    xtarget     = loc_vec[0][Nx]
    ytarget     = loc_vec[1][Ny]
    ztarget     = loc_vec[2][Nz]
    #---------------------------------------------------------------------#
    # time loop                                                           #
    #---------------------------------------------------------------------#
    pcount = 0
    for count, t in enumerate(range(end, start, -1)):
        dt  = time[t] - time[t-1]
        for k in range(0,64):
            for j in range(0,64):
                for i in range(0,64):
                    #-----------------------------------------------------#
                    # finding new points and distance                     #
                    #-----------------------------------------------------#
                    xnew        = loc_vec[0][i] + vel[0][i,j,k,t-1]*dt
                    ynew        = loc_vec[1][j] + vel[1][i,j,k,t-1]*dt
                    znew        = loc_vec[2][k] + vel[2][i,j,k,t-1]*dt
                    dist[i,j,k] = distance(xtarget, xnew, ytarget, ynew, ztarget, znew)
        #-----------------------------------------------------------------#
        # setting new targets and storing points                          #
        #----------------------------------------------------------------#
        [val, Nx, Ny, Nz]   = find_min3D(dist)
        xpt[count]  = Nx
        ypt[count]  = Ny
        zpt[count]  = Nz
        xtarget     = loc_vec[0][Nx]
        ytarget     = loc_vec[1][Ny]
        ztarget     = loc_vec[2][Nz]
        #-----------------------------------------------------------------#
        # print statement                                                 #
        #-----------------------------------------------------------------#
        if pcount > 20:
            print('time step --> %i'        %(count))
            pcount  = 0
        pcount += 1

    return xpt, ypt, zpt


#-------------------------------------------------------------------------#
# Periodic boundary conditions                                            #
#-------------------------------------------------------------------------#
def location_adjustment(
        loc):

    """ Subroutine to account for the periodic boundary conditions """
    #---------------------------------------------------------------------#
    # Domain variables                                                    #
    #---------------------------------------------------------------------#
    PI  = np.pi
    #---------------------------------------------------------------------#
    # Location >  2pi                                                     #
    #---------------------------------------------------------------------#
    if loc > 2.0*PI:
        Loc = loc - 2.0*PI
    #---------------------------------------------------------------------#
    # Location <  0                                                       #
    #---------------------------------------------------------------------#
    elif loc < 0.0:
        Loc = loc + 2.0*PI
    #---------------------------------------------------------------------#
    # 0.0 <= Location <= 2.0*pi                                           #
    #---------------------------------------------------------------------#
    else:
        Loc = loc

    return Loc
#-------------------------------------------------------------------------#
# Subroutine to find specific value                                       #
#-------------------------------------------------------------------------#
def find_value(
        xvec,                   # x vector
        xval):                  # x value

    """ subroutine to find the coordinate value """
    #---------------------------------------------------------------------#
    # Finding each distance between each point                            #
    #---------------------------------------------------------------------#
    dist     = np.zeros(len(xvec))
    for Count, Val in enumerate(xvec):
        dist[Count]  = abs(Val-xval)
    #---------------------------------------------------------------------#
    # Finding the minimum distance                                        #
    #---------------------------------------------------------------------#
    index   = np.unravel_index(np.argmin(dist, axis=None), dist.shape)[0]
    index   = int(index)

    return index
#-------------------------------------------------------------------------#
# Tracking subroutine                                                     #
#-------------------------------------------------------------------------#
def tracking(
        start,                      # start time point
        end,                        # end time point
        Time,                       # time vector
        loc_vec,                    # x, y, and z vectors
        crd_vec,                    # x, y, and z starting coordinates
        vel,                        # vel_x, vel_y, vel_z
        box_size,                   # box size for averaging
        BA_flag = True,             # box average flag
        debug   = False):           # debug flag

    """ Subroutine to track a particular point """
    #---------------------------------------------------------------------#
    # Setting coordinates                                                 #
    #---------------------------------------------------------------------#
    Nx      = crd_vec[0]
    Ny      = crd_vec[1]
    Nz      = crd_vec[2]
    #---------------------------------------------------------------------#
    # Initial conditions                                                  #
    #---------------------------------------------------------------------#
    xt  = loc_vec[0][Nx]
    yt  = loc_vec[1][Ny]
    zt  = loc_vec[2][Nz]
    #---------------------------------------------------------------------#
    # Preallocating time loop variables                                   #
    #---------------------------------------------------------------------#
    x_pts   = np.zeros((end-start)+1)
    y_pts   = np.zeros((end-start)+1)
    z_pts   = np.zeros((end-start)+1)
    #---------------------------------------------------------------------#
    # Time loop                                                           #
    #---------------------------------------------------------------------#
    for t in range(end, start-1, -1):
        #-----------------------------------------------------------------#
        # Storing coordinate values                                       #
        #-----------------------------------------------------------------#
        index           = t-start
        print(index)
        x_pts[index]    = Nx
        y_pts[index]    = Ny
        z_pts[index]    = Nz
        #-----------------------------------------------------------------#
        # previous x, y, and z values                                     #
        #-----------------------------------------------------------------#
        xt_old  = xt
        yt_old  = yt
        zt_old  = zt
        #-----------------------------------------------------------------#
        # Previous Nx, Ny, Nz                                             #
        #-----------------------------------------------------------------#
        Nx_old  = Nx
        Ny_old  = Ny
        Nz_old  = Nz
        #-----------------------------------------------------------------#
        # Finding x, y, and z velocities                                  #
        #-----------------------------------------------------------------#
        if BA_flag is True:
            vel_x   = box_average(vel[0][:,:,:,t],Nx, Ny, Nz, box_size)  # corresponding x-vel.
            vel_y   = box_average(vel[1][:,:,:,t],Nx, Ny, Nz, box_size)  # corresponding x-vel.
            vel_z   = box_average(vel[2][:,:,:,t],Nx, Ny, Nz, box_size)  # corresponding x-vel.
        else:
            vel_x   = star_average(vel[0][:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
            vel_y   = star_average(vel[1][:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
            vel_z   = star_average(vel[2][:,:,:,t],Nx, Ny, Nz)  # corresponding x-vel.
        #-----------------------------------------------------------------#
        # Calculating new locations                                       #
        #-----------------------------------------------------------------#
        dt      = Time[t] - Time[t-1]       # Delta t
        xt      = xt + vel_x*dt             # x location at t^{n-1}
        yt      = yt + vel_y*dt             # y location at t^{n-1}
        zt      = zt + vel_z*dt             # z location at t^{n-1}
        #-----------------------------------------------------------------#
        # Applying periodic boundary conditions                           #
        #-----------------------------------------------------------------#
        xt      = location_adjustment(xt)
        yt      = location_adjustment(yt)
        zt      = location_adjustment(zt)
        #-----------------------------------------------------------------#
        # Finding the x, y, and z coordinate                              #
        #-----------------------------------------------------------------#
        Nx      = find_value(loc_vec[0], xt)
        Ny      = find_value(loc_vec[1], yt)
        Nz      = find_value(loc_vec[2], zt)
        #-----------------------------------------------------------------#
        # Debug print                                                     #
        #-----------------------------------------------------------------#
        if debug is True:
            print('x-values:')
            print('\t\t\ttime=%.5f\t\t\tdt=%.5e\t\t\tvel_{x}=%.5f\t\t\tx_{n+1}=%.5f\t\t\tx_{n}=%.5f\t\t\tNx_{n+1}=%i\t\t\tNx_{n}=%i'\
                    %(t,dt,vel_x,xt_old,xt,Nx_old,Nx))
            print('y-values:')
            print('\t\t\ttime=%.5f\t\t\tdt=%.5e\t\t\tvel_{y}=%.5f\t\t\ty_{n+1}=%.5f\t\t\ty_{n}=%.5f\t\t\tNy_{n+1}=%i\t\t\tNy_{n}=%i'\
                    %(t,dt,vel_y,yt_old,yt,Ny_old,Ny))
            print('z-values:')
            print('\t\t\ttime=%.5f\t\t\tdt=%.5e\t\t\tvel_{z}=%.5f\t\t\tz_{n+1}=%.5f\t\t\tz_{n}=%.5f\t\t\tNz_{n+1}=%i\t\t\tNz_{n}=%i'\
                    %(t,dt,vel_z,zt_old,zt,Nz_old,Nz))
            print('\n\n')

    return x_pts, y_pts, z_pts
#-------------------------------------------------------------------------#
# Star averaging                                                          #
#-------------------------------------------------------------------------#
def star_average(
        field,
        Nx,
        Ny,
        Nz):

    """ Calculating a star average of a single point in the field """
    #---------------------------------------------------------------------#
    # Finding the values of seven points                                  #
    #---------------------------------------------------------------------#
    Vals    = np.zeros(7)
    Vals[0] = field[Nx,Ny,Nz]
    Vals[1] = field[coordinate_adjustment(Nx+1),Ny,Nz]
    Vals[2] = field[coordinate_adjustment(Nx-1),Ny,Nz]
    Vals[3] = field[Nx,coordinate_adjustment(Ny+1),Nz]
    Vals[4] = field[Nx,coordinate_adjustment(Ny-1),Nz]
    Vals[5] = field[Nx,Ny,coordinate_adjustment(Nz+1)]
    Vals[6] = field[Nx,Ny,coordinate_adjustment(Nz-1)]
    Value   = np.mean(Vals)

    return Value
#-------------------------------------------------------------------------#
# Periodic boundary conditions for coordinates                            #
#-------------------------------------------------------------------------#
def coordinate_adjustment(
        coor):

    """ Adjusting the coordinate to range from 0 <= coor <= 63 """
    #---------------------------------------------------------------------#
    # coordinates > 63                                                    #
    #---------------------------------------------------------------------#
    if coor > 63:
        coor = coor - 63
    elif coor < 0:
        coor = 63 + coor
    else:
        coor = coor

    return coor
#-------------------------------------------------------------------------#
# Box average                                                             #
#-------------------------------------------------------------------------#
def box_average(
        field,
        Nx,
        Ny,
        Nz,
        num,
        debug   = False):

    """ Calculating the box average using a different methods """
    #---------------------------------------------------------------------#
    # Setting the box limits                                              #
    #---------------------------------------------------------------------#
    num = int((num-1)/2)
    #---------------------------------------------------------------------#
    # Looping over the box points                                         #
    #---------------------------------------------------------------------#
    Vals    = 0.0
    count   = 0
    for K in range(-num,num+1):
        for J in range(-num, num+1):
            for I in range(-num,num+1):
                Vals    += field[\
                                coordinate_adjustment(Nx+I),\
                                coordinate_adjustment(Ny+J),\
                                coordinate_adjustment(Nz+K)]
                count   += 1
    #---------------------------------------------------------------------#
    # Debug print                                                         #
    #---------------------------------------------------------------------#
    if debug is True:
        print('c=%.2f'      %(count))
        print('Vals=%.5f'   %(Vals))
    #---------------------------------------------------------------------#
    # Taking the average                                                  #
    #---------------------------------------------------------------------#
    Vals    = Vals/count

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
    DS  = True
    print('Loading data:')
    if DS is True:
        time    = np.load(data_path + 'time-2.npy')
        print('\ttime')
        u_x     = np.load(data_path + 'velocity1-2.npy')
        print('\tx-veloccity')
        u_y     = np.load(data_path + 'velocity2-2.npy')
        print('\ty-velocity')
        u_z     = np.load(data_path + 'velocity3-2.npy')
        print('\tz-velocity')
    else:
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
    # Defining final time                                                 #
    #---------------------------------------------------------------------#
    tfinal          = 99
    #---------------------------------------------------------------------#
    # Generating contours plots starting at T=t_{f}                       #
    #---------------------------------------------------------------------#
    if DS is not True:
        [val, X, Y, Z]  = find_max3D(Dk[:,:,:,tfinal])
        for i in range(tfinal, tfinal-6, -1):
            val = find_max3D(Dk[:,:,:,i])[0]
            print('maximum value=%.7f\tx-loc=%.5f\ty-loc=%.5f\tz-loc=%.5f'\
                        %(val, X, Y, Z))
            #contour_gen(Dk, X, Y, Z, i)
    #---------------------------------------------------------------------#
    # Tracking                                                            #
    #---------------------------------------------------------------------#
    [X,Y,Z] = tracking2(tfinal-90, tfinal, time, (x,y,z), [12,60,23], [u_x, u_y,\
                        u_z])
    #---------------------------------------------------------------------#
    # Writing data                                                        #
    #---------------------------------------------------------------------#
    string = ''
    for i, nx in enumerate(X):
        string  += '%i\t%i\t%i\n'               %(X[i], Y[i], Z[i])
    f   = open(data_path + 'box-average-7-pts.txt',  'w')
    f.write(string)
    f.close()
    #---------------------------------------------------------------------#
    # Example plot                                                        #
    #---------------------------------------------------------------------#
    plt.plot(time, u_x[31,31,31,:], 'r', lw=1.5)
    plt.grid(True)
    plt.show()

    print('**** Successful run ****')
    sys.exit(0)

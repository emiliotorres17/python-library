#!/usr/bin/env python3
#-------------------------------------------------------------------------#
# Python packages                                                         #
#-------------------------------------------------------------------------#
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,\
                                AutoMinorLocator
#-------------------------------------------------------------------------#
# Font settings                                                           #
#-------------------------------------------------------------------------#
def plot_setting(
        major_tick  = False,
        minor_tick  = False,
        tick_format = False,
        y_tick      = False):
    #---------------------------------------------------------------------#
    # Font sizes                                                          #
    #---------------------------------------------------------------------#
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    SMALL_SIZE = 14
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 12
    plt.rc('font',      size=SMALL_SIZE)            # controls default text sizes
    plt.rc('axes',      titlesize=SMALL_SIZE)       # fontsize of the axes title
    plt.rc('axes',      labelsize=MEDIUM_SIZE)      # fontsize of the x and y labels
    plt.rc('xtick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    plt.rc('ytick',     labelsize=SMALL_SIZE)       # fontsize of the tick labels
    plt.rc('legend',    fontsize=SMALL_SIZE)        # legend fontsize
    plt.rc('figure',    titlesize=BIGGER_SIZE)      # fontsize of the figure title
    #---------------------------------------------------------------------#
    # Setting tick marks                                                  #
    #---------------------------------------------------------------------#
    #fig, ax = plt.subplots()
    #if major_tick is not False:
    #    ax.xaxis.set_major_locator(MultipleLocator(major_tick))
    #    ax.xaxis.set_major_formatter(FormatStrFormatter(tick_format))
    #    ax.xaxis.set_minor_locator(MultipleLocator(minor_tick))
    #if y_tick is not False:
    #    ax.yaxis.set_major_locator(MultipleLocator(major_tick))
    #    ax.yaxis.set_major_formatter(FormatStrFormatter(tick_format))
    #    ax.yaxis.set_minor_locator(MultipleLocator(minor_tick))
    return
#-------------------------------------------------------------------------#
# Miscellaneous                                                           #
#-------------------------------------------------------------------------#
#def misc(
#    #---------------------------------------------------------------------#
#    # Tick mark settings                                                  #
#    #---------------------------------------------------------------------#
#    fig, ax = plt.subplots()
#    if major_tick is not False:
#        ax.xaxis.set_major_locator(MultipleLocator(major_tick))
#        ax.xaxis.set_major_formatter(FormatStrFormatter(tick_format))
#        ax.xaxis.set_minor_locator(MultipleLocator(minor_tick))
#    if y_tick is not False:
#        ax.yaxis.set_major_locator(MultipleLocator(major_tick))
#        ax.yaxis.set_major_formatter(FormatStrFormatter(tick_format))
#        ax.yaxis.set_minor_locator(MultipleLocator(minor_tick))
#    #---------------------------------------------------------------------#
#    # Plot example for plotting comparison                                #
#    #---------------------------------------------------------------------#
#    fig, ax1    = plt.subplots(1,1)
#    ax1.plot(tvals, error1, 'r', lw=1.5, label='$x_{1}$')
#    ax1.plot(tvals, error2, 'b', lw=1.5, label='$x_{2}$')
#    ax1.plot(tvals, error3, 'g', lw=1.5, label='$x_{3}$')
#    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,1))
#    plt.grid(True)
#    plt.ylabel('Error')
#    plt.xlabel('Time step')
#    plt.legend(loc=0)
#    ax1.invert_xaxis()
#    plt.savefig(media_path + 'error-64.png')
#    plt.close()
if __name__ == '__main__':
    print('Has no been test')
    sys.exit(0)

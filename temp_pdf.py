#!/usr/bin/env python3
"""=============================================================================
 Purpose:

    **** Purpose of the code goes here ****

 Author:
    Emilio E. Torres
============================================================================="""
#==============================================================================#
# Preamble                                                                     #
#==============================================================================#
#------------------------------------------------------------------------------#
# Python packages                                                              #
#------------------------------------------------------------------------------#
import sys
import os
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------#
# User packages                                                                #
#------------------------------------------------------------------------------#

# user packages go here

#==============================================================================#
# User defined functions                                                       #
#==============================================================================#

# user defined functions go here

#==============================================================================#
# Main                                                                         #
#==============================================================================#
if __name__ == '__main__':
    #--------------------------------------------------------------------------#
    # Main preamble                                                            #
    #--------------------------------------------------------------------------#
    call(['clear'])
    sep         = os.sep
    pwd         = os.getcwd()
    data_path   = pwd + '%cdata%c'              %(sep, sep)
    media_path  = pwd + '%cmedia%c'             %(sep, sep)
    #--------------------------------------------------------------------------#
    # Main execution                                                           #
    #--------------------------------------------------------------------------#
    tmp         = [15.841, 16.612, 17.568, 15.812, 17.458, 16.846, 15.855,\
                    17.097, 15.767, 16.906, 16.550]
    values      = np.array(tmp, dtype=float)
    n           = values.size
    mu          = values.mean()
    sigma       = values.std(ddof=1)
    xmin, xmax  = values.min(), values.max()
    xpad        = 0.2*(xmax-xmin) if xmax > xmin else 0.5
    x           = np.linspace(xmin-xpad, xmax+xpad, 800)
    pdf         = (1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2.0)
    bins        = max(5, int(round(np.sqrt(n)*2)))
    fig, ax     = plt.subplots()
    ax.hist(values, bins=bins, density=True, alpha=0.5, edgecolor='black',)
    ax.plot(x, pdf, linewidth=2)
    
    edges = np.histogram_bin_edges(values, bins=bins)

    # bin index 0..bins-1 for each point
    idx = np.digitize(values, edges, right=False) - 1
    idx = np.clip(idx, 0, bins-1)  # keep max in last bin
    
    print(f"bins = {bins}")
    print("edges =", edges)
    
    for i in range(bins):
        left, right = edges[i], edges[i+1]
        bracket_r = "]" if i == bins-1 else ")"
        # +1 so you see point numbers 1..n (instead of 0..n-1)
        point_nums = (np.where(idx == i)[0]).tolist()
        print(f"{i+1}: [{left:.6f}, {right:.6f}{bracket_r}  points={point_nums}  count={len(point_nums)}")

    plt.show()
    


    # main script goes here with

    print('\n **** Successful run ****\n')
    sys.exit(0)

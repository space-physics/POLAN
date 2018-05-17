#!/usr/bin/env python
"""
This example is just a start.
"""
from pathlib import Path
import xarray
import numpy as np
from matplotlib.pyplot import figure,show
#
import runpolan

N = 399
debug = 0
Lmin = 6 # arbitrary, to discard spurious results
infn = Path('examples/in.dat')


def gopolan(infn) -> xarray.Dataset:

    fv = np.zeros(N,dtype=np.float32)
    ht = np.zeros(N,dtype=np.float32)



    fv,ht,qq,fh,dip,start,amode,valley = runpolan.polrun(infn,N,debug)

    height = []
    freqMHz = []
    fr = []; he = []
    for i,(f,h) in enumerate(zip(fv,ht)):
        if len(freqMHz)==0 and f == 0.: #find first value
            continue
        elif f==0. or h==0. or f<fv[i-1] or h<ht[i-1] or f<1. or f>30:
            if len(fr) >= Lmin:
                freqMHz.append(fr)
                height.append(he)
            fr=[]; he=[]
        else:
            fr.append(f); he.append(h)


    iono = {'fv':freqMHz,
            'height':height,
            'dip':dip}

    return iono



def plotiono(iono):

    ax = figure().gca()

    for i,(fv,ht) in enumerate(zip(iono['fv'], iono['height'])):
        ax.plot(fv, ht,label=f'trace {i}')

    ax.set_xlabel('frequency [MHz]')
    ax.set_ylabel('height [km]')
    ax.set_title(f'magnetic dip angle [deg.] {iono["dip"]}')

    ax.grid()

    ax.legend()

if __name__ == '__main__':

    iono = gopolan(infn)


    plotiono(iono)

    show()

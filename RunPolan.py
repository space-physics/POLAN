#!/usr/bin/env python
from pathlib import Path
import numpy as np
#
import runpolan

N = 399
debug = 0
valley= 0.
amode = 0.
start = -1.  #-1 uses a direct start, from the first scaled point.
dip = 20.
fh = -1.

fv = np.zeros(N,dtype=np.float32)
ht = np.zeros(N,dtype=np.float32)

infn = Path('examples/in.asc')

dat = np.loadtxt(infn).reshape(-1,2)
Nr = dat.shape[0]


fv[:Nr] = dat[:,0]
ht[:Nr] = dat[:,1]


qq = runpolan.polan(fv=fv, ht=ht, fb=fh, dip=dip, start=start, amode=amode, valley=valley, list=debug)

print('magnetic dip angle [degrees]',dip)
print(fv)

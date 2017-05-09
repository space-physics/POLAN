#!/usr/bin/env python
from pathlib import Path
import numpy as np
#
import runpolan

N = 399
debug = 0
valley= 0
amode = 0
start = -1.
dip = 20.
fh = -1.

fv = np.zeros(N,dtype=np.float32)
ht = np.zeros(N,dtype=np.float32)

infn = Path('examples/in.asc')

dat = np.loadtxt(infn).reshape(-1,2)
Nr = dat.shape[0]


fv[:Nr] = dat[:,0]
ht[:Nr] = dat[:,1]


qq = runpolan.polan(fv,ht,fh, dip, start, amode, valley, debug)

from typing import List, Dict
import numpy as np
from pathlib import Path

import runpolan

N = 399
debug = 0
Lmin = 6  # arbitrary, to discard spurious results


def gopolan(infn: Path) -> Dict[str, List[float]]:

    fv = np.zeros(N, dtype=np.float32)
    ht = np.zeros(N, dtype=np.float32)

    fv, ht, qq, fh, dip, start, amode, valley = runpolan.polrun(infn, N, debug)

    height = []
    freqMHz: List[float] = []
    fr: List[float] = []
    he: List[float] = []
    for i, (f, h) in enumerate(zip(fv, ht)):
        if len(freqMHz) == 0 and f == 0.:  # find first value
            continue
        elif f == 0. or h == 0. or f < fv[i-1] or h < ht[i-1] or f < 1. or f > 30:
            if len(fr) >= Lmin:
                freqMHz.append(fr)  # type: ignore
                height.append(he)
            fr = []
            he = []
        else:
            fr.append(f)
            he.append(h)

    iono = {'fv': freqMHz,
            'height': height,
            'dip': dip}

    return iono

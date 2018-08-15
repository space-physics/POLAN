#!/usr/bin/env python
"""
This example is just a start.
"""

from pathlib import Path
from matplotlib.pyplot import show

import polan

infn = Path('examples/in.dat')


if __name__ == '__main__':

    iono = polan.gopolan(infn)

    polan.plotiono(iono)

    show()

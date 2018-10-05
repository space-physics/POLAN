#!/usr/bin/env python
"""
This example is just a start.
"""

from pathlib import Path
from matplotlib.pyplot import show

import polan
import polan.plots as pp

infn = Path('examples/in.dat')


def main():
    iono = polan.gopolan(infn)

    pp.plotiono(iono)

    show()


if __name__ == '__main__':
    main()

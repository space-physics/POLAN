import polan
import polan.plot

from pathlib import Path


infn = Path(__file__).parent / "in.dat"

iono = polan.gopolan(infn)
fig = polan.plot.altitude_plot(iono)
fig.show()

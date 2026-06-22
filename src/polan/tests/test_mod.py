from pathlib import Path
import polan
from pytest import approx

R = Path(__file__).resolve().parents[1]
infn = R / "examples/in.dat"


def test_basic():

    iono = polan.gopolan(infn)

    assert iono["dip"] == 0.0
    assert iono["fv"][0][0] == approx(1.774)

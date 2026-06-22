from pathlib import Path
import polan
import pytest
from pytest import approx

R = Path(__file__).resolve().parents[3]
infn = R / "examples/in.dat"


def test_basic():

    iono = polan.gopolan(infn)

    assert iono["dip"] == 0.0
    assert iono["fv"][0][0] == approx(1.774)


def test_missing_input_raises_file_not_found():

    missing = R / "examples/does-not-exist.dat"
    with pytest.raises(FileNotFoundError):
        polan.gopolan(missing)

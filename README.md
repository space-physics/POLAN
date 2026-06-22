# POLAN

[![ci](https://github.com/space-physics/POLAN/actions/workflows/ci.yml/badge.svg)](https://github.com/space-physics/POLAN/actions/workflows/ci.yml)

POLAN is a classic Fortran program by J. E. Titheridge used to calculate real-height
profiles from chirp ionosonde data from the ionosphere.
The
[POLAN code](https://web.archive.org/web/20231111062759/https://sws.bom.gov.au/IPSHosted/INAG/uag_93/uag_93.html)
was updated to compile on modern PCs.

We provide a Python interface to Fortran POLAN that is installed as usual with pip:

```sh
git clone https://github.com/space-physics/POLAN.git

pip install ./POLAN
```

Optionally one can install POLAN Python interface in
[editable module mode](https://pip.pypa.io/en/stable/topics/local-project-installs/#editable-installs),
which requires disabling build isolation
[disabled](https://pip.pypa.io/en/stable/reference/build-system/#disabling-build-isolation):

```sh
pip install numpy meson meson-python ninja

pip install -e --no-build-isolation ./POLAN
```

## CMake build (optional)

For those who wish to use POLAN as a standalone Fortran program, it can be built using CMake.

Build POLAN:

```sh
cmake --workflow default
```

Run an example POLAN calculation

```sh
ctest --test-dir build -V
```

this creates a big output text file `out.dat`

See [Readme_polan.md](./Readme_polan.md) for more details on the POLAN code and how to run it.

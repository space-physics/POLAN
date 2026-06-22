py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Discover NumPy include directory at configure time
numpy_include = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

inc_numpy = include_directories(numpy_include)

fortran_sources = files(
  'polan.f',
  'polmis.f',
  'polsin.f',
  'polsub.f',
  'polrunsub.f',
)

# Generate C API wrapper using numpy.f2py and move outputs to build directory
f2py_wrapper = custom_target(
  'runpolan-f2py-wrapper',
  input: fortran_sources,
  output: ['runpolanmodule.c', 'runpolan-f2pywrappers.f'],
  command: [
    py, '-c',
    'import pathlib, shutil, subprocess, sys; outdir = pathlib.Path(sys.argv[1]); subprocess.run([sys.executable, "-m", "numpy.f2py", *sys.argv[2:], "-m", "runpolan", "--lower"], check=True); shutil.move("runpolanmodule.c", outdir / "runpolanmodule.c"); shutil.move("runpolan-f2pywrappers.f", outdir / "runpolan-f2pywrappers.f")',
    '@OUTDIR@',
    '@INPUT@',
  ],
)

# Copy f2py runtime support headers from NumPy to stable build location
f2py_runtime = custom_target(
  'runpolan-f2py-runtime',
  output: ['fortranobject.c', 'fortranobject.h'],
  command: [
    py, '-c',
    'import shutil, sys, pathlib; import numpy.f2py as f2py; src = pathlib.Path(f2py.get_include()); out = pathlib.Path(sys.argv[1]); shutil.copyfile(src / "fortranobject.c", out / "fortranobject.c"); shutil.copyfile(src / "fortranobject.h", out / "fortranobject.h")',
    '@OUTDIR@',
  ],
)

# Build extension module
py.extension_module(
  'runpolan',
  [
    fortran_sources,
    f2py_wrapper,
    f2py_runtime[0],
  ],
  dependencies: py_dep,
  include_directories: [inc_numpy],
  link_language: 'fortran',
  install: true,
  subdir: 'polan',
)

py.install_sources('polan/__init__.py', subdir: 'polan')

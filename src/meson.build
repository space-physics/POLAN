py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# NumPy and f2py include directories (canonical f2py + Meson pattern)
# https://numpy.org/doc/stable/f2py/buildtools/meson.html
incdir_numpy = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true
).stdout().strip()

incdir_f2py = run_command(
  py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

fortran_sources = files(
  'polan.f',
  'polmis.f',
  'polsin.f',
  'polsub.f',
  'polrunsub.f',
)

# Generate the C API wrapper from Fortran sources with numpy.f2py.
# '--build-dir @OUTDIR@' makes f2py write its outputs where Meson expects them.
runpolan_source = custom_target(
  'runpolanmodule.c',
  input: fortran_sources,
  output: ['runpolanmodule.c', 'runpolan-f2pywrappers.f'],
  command: [
    py, '-m', 'numpy.f2py',
    '@INPUT@',
    '-m', 'runpolan',
    '--lower',
    '--build-dir', '@OUTDIR@',
  ],
)

# On Windows with gfortran (MinGW), statically link the Fortran runtime so the
# extension module is self-contained and does not need the MinGW runtime DLLs
# (libgfortran, libgcc, libquadmath, libwinpthread) on PATH at import time.
fortran_compiler = meson.get_compiler('fortran')
extra_link_args = []
if host_machine.system() == 'windows' and fortran_compiler.get_id() == 'gcc'
  extra_link_args += [
    '-static-libgcc',
    '-static-libgfortran',
    '-static-libquadmath',
    '-Wl,-Bstatic,--whole-archive',
    '-lwinpthread',
    '-Wl,--no-whole-archive',
  ]
endif

# Build the extension module. fortranobject.c is referenced directly from the
# installed NumPy f2py include directory; no copying required.
py.extension_module(
  'runpolan',
  [
    fortran_sources,
    runpolan_source,
    incdir_f2py / 'fortranobject.c',
  ],
  include_directories: inc_np,
  dependencies: py_dep,
  link_language: 'fortran',
  link_args: extra_link_args,
  install: true,
  subdir: 'polan',
)

py.install_sources('polan/__init__.py', subdir: 'polan')

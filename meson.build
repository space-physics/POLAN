project(
  'polan',
  ['c', 'fortran'],
  version: '1.1.0',
  meson_version: '>=1.2.0'
)

# Enable legacy Fortran standard for GFortran compatibility with fixed-form code
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_project_arguments('-std=legacy', language: 'fortran')
endif

subdir('src')

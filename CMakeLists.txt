cmake_minimum_required(VERSION 2.8.12)
project( polan Fortran )

add_compile_options(-mtune=native -Wno-line-truncation)

add_executable(polan polrun.f polan.f polmis.f polsin.f polsub.f)

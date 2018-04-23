cmake_minimum_required(VERSION 3.10)
project(polan Fortran)
enable_testing()

add_compile_options(-mtune=native -g -O3) 

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_compile_options(-Wno-line-truncation)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  link_libraries(-static-flang-libs)
endif()

add_executable(polan polrun.f polan.f polmis.f polsin.f polsub.f)

add_test(NAME basic COMMAND ./polan ${CMAKE_CURRENT_SOURCE_DIR}/examples/in.dat)

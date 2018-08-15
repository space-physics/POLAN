cmake_minimum_required(VERSION 3.3)
project(polan Fortran)
enable_testing()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_compile_options(-mtune=native -Werror=array-bounds -Wno-line-truncation)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
  add_compile_options(-Mallocatable=03)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang) 
  add_compile_options(-Mallocatable=03)
  link_libraries(-static-flang-libs)
endif()

add_executable(polan src/polrun.f src/polan.f src/polmis.f src/polsin.f src/polsub.f)

add_test(NAME basic COMMAND ./polan ${CMAKE_CURRENT_SOURCE_DIR}/examples/in.dat)

if(NOT WIN32)
  add_test(NAME CompareOut COMMAND bash -c "diff <(tail -n +2 out.dat) <(tail -n +2 ../examples/out.dat)")
  set_tests_properties(CompareOut PROPERTIES FIXTURES_SETUP basic)
endif()

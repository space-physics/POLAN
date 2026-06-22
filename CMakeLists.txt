cmake_minimum_required(VERSION 3.19...4.4)

project(polan LANGUAGES Fortran)

enable_testing()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-std=legacy)
endif()

add_executable(polan)
add_test(NAME basic COMMAND polan ${polan_SOURCE_DIR}/examples/in.dat)

add_subdirectory(src)

if(UNIX)
  add_test(NAME CompareOut COMMAND bash -c "diff <(tail -n +2 out.dat) <(tail -n +2 ../examples/out.dat)")
  set_tests_properties(CompareOut PROPERTIES FIXTURES_SETUP basic)
endif()

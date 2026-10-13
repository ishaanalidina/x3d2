# Version 4.2 is required at minimum when using the iterative solver. This version fixes builds failing due to increased circular dependencies in later PETSc versions (>3.23). The circular dependency issues have been highlighted in PETSc issues #1738 and #1771, while the CMake update version requirement has been seen in CMake issue #27238.

cmake_minimum_required(VERSION 4.2)
project(x3d2 LANGUAGES Fortran)
enable_testing()

set(CMAKE_DEPENDS_USE_COMPILER TRUE)

#
# Set the Poisson solver choice
#
set(POISSON_SOLVER FFT CACHE STRING
  "Select the Poisson solver: FFT or ITER")
set(WITH_2DECOMPFFT ON CACHE BOOL
  "Enable Poisson based FFT solver on the OpenMP backend.")

add_subdirectory(src)
add_subdirectory(tests)

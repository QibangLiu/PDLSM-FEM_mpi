# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIRS - where to find mkl.h, etc.
#  MKL_LIBRARIES    - List of libraries when using mkl.
#  MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl_dfti.h HINTS $ENV{MKLROOT}/include)
#find_path(FFTW_INCLUDE_DIR fftw.h HINTS $ENV{MKLROOT}/include/fftw)
set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR})
#message(STAUS "==${MKL_INCLUDE_DIR}")

if(WIN32)
    file(GLOB  MKL_LIBRARY 
    #$ENV{MKLROOT}/lib/intel64/mkl_rt.lib
    $ENV{MKLROOT}/lib/intel64/mkl_core.lib
    #$ENV{MKLROOT}/lib/intel64/mkl_intel_ilp64.lib
    $ENV{MKLROOT}/lib/intel64/mkl_intel_lp64.lib
    #$ENV{MKLROOT}/lib/intel64/mkl_intel_thread.lib
    $ENV{MKLROOT}/lib/intel64/mkl_sequential.lib)
elseif(UNIX)
    file(GLOB  MKL_LIBRARY 
    $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.so
    $ENV{MKLROOT}/lib/intel64/libmkl_core.so
    $ENV{MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.so
    $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.so
    $ENV{MKLROOT}/lib/intel64/libmkl_rt.so)
   #find_library(MKL_LIBRARY NAMES mkl_rt  HINTS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
endif()

set(MKL_LIBRARIES ${MKL_LIBRARY})
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARY MKL_INCLUDE_DIR)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY )


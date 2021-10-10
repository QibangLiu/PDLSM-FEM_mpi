# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIRS - where to find mkl.h, etc.
#  MKL_LIBRARIES    - List of libraries when using mkl.
#  MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl.h HINTS $ENV{MKLROOT}/include)
#find_path(FFTW_INCLUDE_DIR fftw.h HINTS $ENV{MKLROOT}/include/fftw)
#message(STAUS "==$ENV{MKLROOT}")

if(WIN32)
    file(GLOB  MKL_LIBRARY 
    #$ENV{MKLROOT}/lib/intel64/mkl_rt.lib
    $ENV{MKLROOT}/lib/intel64/mkl_intel_ilp64.lib
    $ENV{MKLROOT}/lib/intel64/mkl_core.lib
    $ENV{MKLROOT}/lib/intel64/mkl_sequential.lib
    $ENV{MKLROOT}/lib/intel64/mkl_blacs_msmpi_ilp64.lib)
elseif(UNIX AND NOT APPLE)
    file(GLOB  MKL_LIBRARY 
    $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.so
    $ENV{MKLROOT}/lib/intel64/libmkl_core.so
    $ENV{MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.so
    $ENV{MKLROOT}/lib/intel64/libmkl_sequential.so)
   #find_library(MKL_LIBRARY NAMES mkl_rt  HINTS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
elseif(APPLE)
    file(GLOB  MKL_LIBRARY 
    $ENV{MKLROOT}/lib/intel64/libmkl_intel_ilp64.dylib
    $ENV{MKLROOT}/lib/intel64/libmkl_sequential.dylib
    $ENV{MKLROOT}/lib/intel64/libmkl_core.dylib
    $ENV{MKLROOT}/lib/intel64/lmkl_blacs_mpich_ilp64.dylib)
endif()

set(MKL_LIBRARIES ${MKL_LIBRARY})
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARY MKL_INCLUDE_DIR)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY )


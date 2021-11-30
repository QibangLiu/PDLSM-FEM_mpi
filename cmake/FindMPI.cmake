# - Find MPI
# Find the native MPI headers and libraries.
#
#  MPI_INCLUDE_DIRS - where to find mpi.h, etc.
#  MPI_LIBRARIES    - List of libraries when using mpi.
#  MPI_FOUND        - True if mpi found.
#

#==========find includ dir===========
if(WIN32)
   find_path(MPI_INCLUDE_DIR mpi.h HINTS $ENV{MPIROOT}/Include $ENV{MPIROOT}/include)
elseif(UNIX)
   find_path(MPI_INCLUDE_DIR mpi.h HINTS $ENV{MPIROOT}/include)
endif()
#============find libraries
if(WIN32)
   find_library(MPI_LIBRARY NAMES msmpi.lib  HINTS $ENV{MPIROOT}/Lib/x64)
elseif(UNIX AND NOT APPLE)
   find_library(MPI_LIBRARY NAMES libmpi.so  HINTS $ENV{MPIROOT}/lib/release $ENV{MPIROOT}/lib)
elseif(APPLE)
   find_library(MPI_LIBRARY NAMES libmpi.dylib  HINTS $ENV{MPIROOT}/lib/release $ENV{MPIROOT}/lib)
endif()
#message(STAUS "==$ENV{MPIROOT}")
set(MPI_LIBRARIES ${MPI_LIBRARY})
set(MPI_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MPI_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MPI DEFAULT_MSG MPI_LIBRARY MPI_INCLUDE_DIR)

mark_as_advanced(MPI_INCLUDE_DIR MPI_LIBRARY )

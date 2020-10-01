########################################
# CMake build system
# This file is part of PDLSM-FEM
# Created by Qibang Liu
cmake_minimum_required(VERSION 3.1)
project(PDLSM-FEM CXX)

#====load macro==================
include(${CMAKE_CURRENT_SOURCE_DIR}/Macros/subdirectoryList.cmake)
#=============================
#====inclusion of my code====================#
get_filename_component(PDFEM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../SourceFiles ABSOLUTE)
SUBDIRLIST(PDFEM_HEADFILE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../HeaderFiles) # all subdirs;
include_directories(${PDFEM_HEADFILE_DIR})
include_directories(${FORMS_DIR})
file(GLOB_RECURSE  PDFEM_SOURCES 
    ${PDFEM_SOURCE_DIR}/*.cpp
     ${PDFEM_HEADFILE_DIR}/*.h
    )
#===Create the source groups for source tree with root at CMAKE_CURRENT_SOURCE_DIR.
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/.. FILES ${PDFEM_SOURCES})

#======inclusion of MKL ===========================
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})
#======inclusion of MPI ===========================
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_DIRS})
#=============ADD the executable===========#
if(WIN32)
    add_executable(${PROJECT_NAME} ${PDFEM_SOURCES})
elseif(UNIX)
    add_executable(${PROJECT_NAME} ${PDFEM_SOURCES})
endif()
#======add the target includes============
target_include_directories(${PROJECT_NAME} PUBLIC ${PDFEM_SOURCE_DIR})
#======linking librarys=====================
target_link_libraries(${PROJECT_NAME} 
-Wl,--start-group
${MKL_LIBRARIES} 
-Wl,--end-group
${MPI_LIBRARIES}
)
#===set complie options=====================
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -lgomp -lpthread -lm -ldl -DMKL_ILP64 -m64")
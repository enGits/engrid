# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                      +
# + This file is part of enGrid.                                         +
# +                                                                      +
# + Copyright 2008-2014 enGits GmbH                                      +
# +                                                                      +
# + enGrid is free software: you can redistribute it and/or modify       +
# + it under the terms of the GNU General Public License as published by +
# + the Free Software Foundation, either version 3 of the License, or    +
# + (at your option) any later version.                                  +
# +                                                                      +
# + enGrid is distributed in the hope that it will be useful,            +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
# + GNU General Public License for more details.                         +
# +                                                                      +
# + You should have received a copy of the GNU General Public License    +
# + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
# +                                                                      +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cmake_minimum_required(VERSION 2.8)

# limit configuration types (must be done before project() statement)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited config" FORCE)

project(ENGRID)

set(CMAKE_VERBOSE_MAKEFILE off)

set(ENGRID_MAJOR_VERSION 2)
set(ENGRID_MINOR_VERSION 0)
set(ENGRID_PATCH_VERSION 0)
set(ENGRID_VERSION ${ENGRID_MAJOR_VERSION}.${ENGRID_MINOR_VERSION}.${ENGRID_PATCH_VERSION})

FIND_PACKAGE(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

FIND_PACKAGE(CGAL REQUIRED)
include(${CGAL_USE_FILE})
#include_directories(${VTK_INCLUDE_DIRS})

# Qt4
find_package(Qt4 COMPONENTS QtCore REQUIRED)
find_package(Qt4 COMPONENTS QtGui REQUIRED)
find_package(Qt4 COMPONENTS QtXml REQUIRED)

#include(${QT_USE_FILE})
include_directories(${QT_QTCORE_INCLUDE_DIR})
include_directories(${QT_QTGUI_INCLUDE_DIR})
include_directories(${QT_QTXML_INCLUDE_DIR})

# set include directory, add src directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tetgen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libengrid)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/libengrid)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_definitions(-DNDEBUG)
  add_definitions(-DQT_NO_DEBUG)
  remove_definitions(-DQT_DEBUG)
else()
  add_definitions(-DQT_DEBUG)
endif()

add_subdirectory(tetgen)
add_subdirectory(libengrid)

add_executable(engrid main.cpp)

# CGAL
#find_path(CGAL_INCLUDE_PATH CGAL/Exact_predicates_inexact_constructions_kernel.h /usr/include /usr/local/include)
#add_definitions(-DUSE_CGAL_PREDICATES)
#include_directories(${CGAL_INCLUDE_PATH})
#find_path (BOOST_INCLUDE_PATH any.hpp /usr/include/boost /usr/local/include/boost)
#include_directories(${BOOST_INCLUDE_PATH})
#find_library(GMP_LIBRARY
#             NAMES libgmp.so
#             PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
#find_library(CGAL_LIBRARY
#             NAMES libCGAL.so
#             PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
#find_library(CGAL_CORE_LIBRARY
#             NAMES libCGAL_Core.so
#             PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
#find_library(BOOST_SYSTEM_LIBRARY
#             NAMES libboost_system.so
#             PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)

target_link_libraries(engrid libengrid)
target_link_libraries(engrid tet)
target_link_libraries(engrid ${QT_LIBRARIES})
target_link_libraries(engrid ${VTK_LIBRARIES})
target_link_libraries(engrid Qt4::QtCore)
target_link_libraries(engrid Qt4::QtGui)
target_link_libraries(engrid Qt4::QtXml)
#target_link_libraries(engrid ${GMP_LIBRARY})
#target_link_libraries(engrid ${BOOST_SYSTEM_LIBRARY})
#target_link_libraries(engrid ${GGAL_CORE_LIBRARY})
#target_link_libraries(engrid ${GGAL_LIBRARY})

install(TARGETS engrid DESTINATION bin)

add_dependencies(engrid tetgen libengrid)



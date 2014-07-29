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

exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "describe"
    OUTPUT_VARIABLE GIT_SHA1)
string(REGEX MATCH "-g.*$" GIT_SHA1 ${GIT_SHA1})
string(REGEX REPLACE "[-g]" "" GIT_SHA1 ${GIT_SHA1} )
add_definitions( -DGIT_SHA1="${GIT_SHA1}" )

exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "symbolic-ref --short HEAD"
    OUTPUT_VARIABLE GIT_BRANCH)
add_definitions( -DGIT_BRANCH="${GIT_BRANCH}" )

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

option(PROFILING "add -g for profiling even in RELEASE configuration?" OFF)

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  add_definitions(-DNDEBUG)
  add_definitions(-DQT_NO_DEBUG)
  remove_definitions(-DQT_DEBUG)
  if (PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  endif()
else()
  add_definitions(-DQT_DEBUG)
endif()

add_subdirectory(tetgen)
add_subdirectory(libengrid)

add_executable(engrid main.cpp)

target_link_libraries(engrid libengrid)
target_link_libraries(engrid tet)
target_link_libraries(engrid ${QT_LIBRARIES})
target_link_libraries(engrid ${VTK_LIBRARIES})
target_link_libraries(engrid Qt4::QtCore)
target_link_libraries(engrid Qt4::QtGui)
target_link_libraries(engrid Qt4::QtXml)

install(TARGETS engrid DESTINATION bin)

install(
  FILES engrid.bash
  PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ
  DESTINATION bin
)

add_dependencies(engrid tetgen libengrid)



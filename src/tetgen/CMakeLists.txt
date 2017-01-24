cmake_minimum_required(VERSION 2.6)

add_definitions(-DNDEBUG)

# Add an executable to the project using the specified source files.
add_executable(tetgen tetgen.cxx predicates.cxx)

# CGAL
find_path(CGAL_INCLUDE_PATH CGAL/Exact_predicates_inexact_constructions_kernel.h /usr/include /usr/local/include)
if(CGAL_INCLUDE_PATH)
  add_definitions(-DUSE_CGAL_PREDICATES)
  find_path (BOOST_INCLUDE_PATH any.hpp /usr/include/boost /usr/local/include/boost /opt/local/include/boost /macports/include)
  include_directories(${BOOST_INCLUDE_PATH})
  #FIND_PACKAGE(Boost)
  #if (BOOST_FOUND)
  #  include_directories(${BOOST_INCLUDE_PATH})
  #endif (BOOST_FOUND)
  find_library(GMP_LIBRARY
               NAMES libgmp.so libgmp.dylib
               PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /opt/local/lib)
  target_link_libraries(tetgen ${GMP_LIBRARY})
endif()

#Add a library to the project using the specified source files. 
# In Linux/Unix, it will creates the libtet.a
add_library(tet STATIC tetgen.cxx predicates.cxx)

#Set properties on a target. 
#We use this here to set -DTETLIBRARY for when compiling the
#library
set_target_properties(tet PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)

install(TARGETS tetgen DESTINATION bin)

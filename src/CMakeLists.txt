#-----------------------------------------------------------------------------
#
# Cleaver2 - A MultiMaterial Tetrahedral Mesher
# -- CMakeList
#
# Primary Author:        Jonathan R. Bronson (bronson@sci.utah.edu)
# Contributing Authors:  Joshua   A. Levine  (levine@sci.utah.edu)
#                        Shankar     Sastry  (sastry@sci.utah.edu)
#
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  Copyright (C) 2011, 2012, 2013, Jonathan Bronson
#  Scientific Computing and Imaging Institute
#  University of Utah
#
#
#  Permission is hereby  granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to  use, copy, modify,  merge, publish, distribute, sublicense,
#  and/or  sell  copies of  the Software, and  to permit persons to whom the
#  Software is  furnished to  do so,  subject  to  the following conditions:
#
#  The above copyright notice  and  this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING  BUT NOT  LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A  PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN  NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER  IN  AN  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT  OF  OR IN  CONNECTION  WITH THE  SOFTWARE  OR  THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#-----------------------------------------------------------------------------

PROJECT(CLEAVER2)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE RelWithDebInfo)
#SET(CMAKE_CXX_FLAGS "-Wall -ansi -Weffc++" )
#SET(CMAKE_DEBUG_FLAGS "-Wall -ansi -Weffc++" )

SET(CMAKE_CXX_FLAGS "-g -std=c++11 ${CMAKE_CXX_FLAGS}" )

if (__APPLE__)
  SET(CMAKE_CXX_FLAGS "--sysroot ${CMAKE_OSX_SYSROOT} ${CMAKE_CXX_FLAGS}")
endif(__APPLE__)

SET(CMAKE_DEBUG_FLAGS "-Wall -ansi -O3" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


# Build Cleaver Library (no dependencies)
ADD_SUBDIRECTORY(lib/teem)
ADD_SUBDIRECTORY(lib/cleaver)
ADD_SUBDIRECTORY(gui)	
ADD_SUBDIRECTORY(cli/mesher)


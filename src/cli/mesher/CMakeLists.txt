PROJECT(cleavercli)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Find QT4, OpenGL
FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(OpenGL)

if(WIN32)
  FIND_PACKAGE(GLEW)
endif()

# Ensure GUI and OpenGL Work
SET( QT_USE_QTGUI TRUE )
SET( QT_USE_QTOPENGL TRUE )
SET( QT_USE_THREAD TRUE )
SET( QT_USE_QTXML TRUE )

# Qt Includes
INCLUDE( ${QT_USE_FILE} )

SET( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/")

#INCLUDE_DIRECTORIES(./Application)

#FILE(GLOB moc_srcs ./Application/*.cpp)
#FILE(GLOB moc_hdrs ./Application/*.h)


# Add library incl1udes
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

# Add Source Files
FILE(GLOB srcs *.cpp TestData/*.cpp)
FILE(GLOB hdrs *.h   TestData/*.h)



#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ./)

#ADD_CUSTOM_TARGET(comp_shaders
#                  DEPENDS ${default_vertex_shader} ${default_fragment_shader})

ADD_EXECUTABLE(cleaver-cli ${srcs} ${hdrs})


TARGET_LINK_LIBRARIES(cleaver-cli teem cleaver )




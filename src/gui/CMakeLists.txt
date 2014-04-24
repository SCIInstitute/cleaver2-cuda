PROJECT(cleavergui)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


# Find QT4, OpenGL
find_package(BLAS QUIET)
FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(OpenGL)

if(WIN32)
  FIND_PACKAGE(GLEW)
endif()

OPTION(BUILD_CLEAVER_APP "Build Cleaver GUI App (REQUIRES TEEM)" ON)

# Ensure GUI and OpenGL Work
SET( QT_USE_QTGUI TRUE )
SET( QT_USE_QTOPENGL TRUE )
SET( QT_USE_THREAD TRUE )
SET( QT_USE_QTXML TRUE )
ADD_DEFINITIONS(-DTETLIBRARY)

# Qt Includes
INCLUDE( ${QT_USE_FILE} )

SET( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/")

INCLUDE_DIRECTORIES(./Application ./Application/DataWidgets)

FILE(GLOB moc_srcs ./Application/*.cpp ./Application/Data/*.cpp ./Application/DataWidgets/*.cpp ./Application/ToolWidgets/*.cpp ./Application/ViewWidgets/*.cpp)
FILE(GLOB moc_hdrs ./Application/*.h   ./Application/Data/*.h   ./Application/DataWidgets/*.h   ./Application/ToolWidgets/*.h   ./Application/ViewWidgets/*.h)

# Add QT Files
FILE(GLOB Forms_UIS ./Application/*.ui ./Application/DataWidgets/*.ui ./Application/ToolWidgets/*.ui ./Application/ViewWidgets/*.ui)
FILE(GLOB RCC_FILES ./Resources/*.qrc)

# Add library incl1udes
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(./Util)

# Add Source Files
FILE(GLOB srcs *.cpp ./Main/*.cpp ./TestData/*.cpp ./Util/*.cpp )
FILE(GLOB hdrs *.h   ./Main/*.h   ./TestData/*.h   ./Util/*.h   )

# Add Shader Files
FILE(GLOB vert ./Application/Shaders/*.vert)
FILE(GLOB frag ./Application/Shaders/*.frag)
FILE(GLOB vf_h ./Application/Shaders/*.h)

#MACRO(compile_shader arg1 arg2)
#    GET_FILENAME_COMPONENT(shader_path ${arg1} ABSOLUTE)
#    GET_FILENAME_COMPONENT(shader_name ${arg1} NAME_WE)
#    GET_FILENAME_COMPONENT(shader_ext ${arg1} EXT)
#    SET(shader_file "${shader_name}${shader_ext}")
#    if(shader_ext STREQUAL ".vert")
#        SET(${arg2} "${shader_name}_vert.h")
#        ADD_CUSTOM_COMMAND(
#            OUTPUT  ${arg2}
#            #MESSAGE(STATUS "Compiling ${arg1} into ${arg2}")
#            COMMAND ${CMAKE_COMMAND} -E copy ${arg1} ${CMAKE_CURRENT_SOURCE_DIR}/${shader_file}
#            COMMAND xxd -i ${shader_file} ${arg2}
#            DEPENDS ${arg1}
#        )
#    elseif(shader_ext STREQUAL ".frag")
#        SET(${arg2} "${shader_name}_frag.h")
#        ADD_CUSTOM_COMMAND(
#            OUTPUT  ${arg2}
#            COMMAND ${CMAKE_COMMAND} -E copy ${arg1} ${CMAKE_CURRENT_SOURCE_DIR}/${shader_file}
#            COMMAND xxd -i ${shader_file} ${arg2}
#            DEPENDS ${arg1}
#        )
#    else()
#        MESSAGE(STATUS "invalid shader given to function Compile_Shader()")
#    endif()
#ENDMACRO()

#COMPILE_SHADER("default.vert" default_vertex_shader)
#MESSAGE(STATUS "default_shader = ${default_vertex_shader}")

#COMPILE_SHADER("default.frag" default_fragment_shader)
#MESSAGE(STATUS "default_shader = ${default_fragment_shader}")


ADD_CUSTOM_COMMAND(
    OUTPUT  default_vert.h
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Application/Shaders/default.vert ${CMAKE_CURRENT_SOURCE_DIR}/default.vert
    COMMAND xxd -i default.vert default_vert.h
    DEPENDS ./Application/Shaders/default.vert
)

ADD_CUSTOM_COMMAND(
    OUTPUT  default_frag.h
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Application/Shaders/default.frag ${CMAKE_CURRENT_SOURCE_DIR}/default.frag
    COMMAND xxd -i default.frag default_frag.h
    DEPENDS ./Application/Shaders/default.frag
)

# Generate QT Stuff
QT4_WRAP_UI(UI_SRC ${Forms_UIS})
QT4_WRAP_CPP(MOC_SRC ${moc_hdrs})
QT4_ADD_RESOURCES(RCC_SRC ${RCC_FILES})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ./)

ADD_CUSTOM_TARGET(comp_shaders
                  DEPENDS ${default_vertex_shader} ${default_fragment_shader})

#include CUDA if we have it
find_package(CUDA)

if(CUDA_FOUND)
	ADD_DEFINITIONS(-DCUDA_FOUND)
	# nvcc flags
	if(CMAKE_COMPILER_IS_GNUCC)
		set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_20,code=sm_20)
	  	string(REPLACE "-std=c++11" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	  	string(REPLACE "-std=c++0x" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	endif()
	cuda_compile(CLEAVER_CUDA_O ./Util/cleaver_cuda.cu ./Util/cleaver_cuda.hh)
	
endif()

SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}" )

if(APPLE)
   ADD_EXECUTABLE(cleaver-gui MACOSX_BUNDLE
       ${srcs} ${hdrs} ${CLEAVER_CUDA_O} ${moc_srcs} ${MOC_SRC} ${UI_SRC} 
       ${RCC_SRC} ${vert} ${frag} ${vf_h})
else()
   ADD_EXECUTABLE(cleaver-gui ${srcs} ${hdrs} ${CLEAVER_CUDA_O} ${moc_srcs}
   ${MOC_SRC} ${UI_SRC} ${RCC_SRC} ${vert} ${frag} ${vf_h})
	   
endif()

if(WIN32)
TARGET_LINK_LIBRARIES(cleaver-gui teem cleaver 
	${OPENGL_LIBRARIES} ${QT_LIBRARIES} ${GLEW_LIBRARIES} ${CUDA_LIBRARIES}) #
else()
TARGET_LINK_LIBRARIES(cleaver-gui teem cleaver 
	${OPENGL_LIBRARIES} ${QT_LIBRARIES} ${CUDA_LIBRARIES}) # stellar
endif()

ADD_DEPENDENCIES(cleaver-gui comp_shaders)



cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME kelvinlets)
project(${PROJECT_NAME})

include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)
include_directories(/Users/xiechulin/code/dependency/eigen-eigen-323c052e1731)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wno-deprecated-declarations")

set(LINK_FLAGS -m64)
add_definitions(-m64)

set(SOURCE_FILES
		KelvinLetsUI.mel
		BrushMPxToolCmd.h
		BrushMPxToolCmd.cpp
		BrushMPxCtx.h
		BrushMPxCtx.cpp
		BrushMPxCtxCmd.h
		BrushMPxCtxCmd.cpp)

# set linking libraries
set(LIBRARIES
		OpenMaya
		OpenMayaUI
		OpenMayaRender
		OpenMayaAnim
		Foundation
		)

build_plugin()




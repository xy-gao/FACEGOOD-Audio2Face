cmake_minimum_required(VERSION 3.10)

project(LPC)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(CMAKE_PREFIX_PATH "${PROJECT_BINARY_DIR}/install/")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

add_compile_definitions(_CRT_SECURE_NO_WARNINGS)

set(FG_EXAMPLE_SRC_LIST)
set(FG_EXAMPLE_INCLUDE_LIST)
set(FG_EXAMPLE_HEADER_LIST)

set(project_inc_dir)
list(APPEND project_inc_dir "${PROJECT_SOURCE_DIR}/src")

foreach(SUBDIR IN LISTS project_inc_dir)
    file(GLOB root_header_files "${SUBDIR}/*.h")
    list(APPEND FG_EXAMPLE_HEADER_LIST "${root_header_files}")    
endforeach()

set(project_src_dir)
list(APPEND project_src_dir "${PROJECT_SOURCE_DIR}/src")

foreach(SUBDIR IN LISTS project_src_dir)
    file(GLOB root_src_files "${SUBDIR}/*.cpp")
    list(APPEND FG_EXAMPLE_SRC_LIST "${root_src_files}")
endforeach()


add_library(LPC SHARED  ${FG_EXAMPLE_SRC_LIST} ${FG_EXAMPLE_HEADER_LIST})

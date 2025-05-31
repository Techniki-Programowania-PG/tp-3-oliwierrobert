cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(sygnaly)
include(FetchContent)

FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)

FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()    
target_compile_features(matplot PUBLIC cxx_std_17)
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(sygnaly main.cpp)
target_link_libraries(sygnaly PUBLIC matplot)






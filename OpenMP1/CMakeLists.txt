cmake_minimum_required(VERSION 3.17)
project(OpenMP1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

set(SOURCE_FILES main.cpp)
add_executable(OpenMP1 ${SOURCE_FILES})
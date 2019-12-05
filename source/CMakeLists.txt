cmake_minimum_required(VERSION 2.8.12)
project(screwIK)

find_package (Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
add_subdirectory(pybind11)
pybind11_add_module(screwIK screwIK.cpp)


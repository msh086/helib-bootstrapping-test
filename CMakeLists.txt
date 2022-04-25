cmake_minimum_required(VERSION 3.21)
project(bootstrapping)

set(CMAKE_CXX_STANDARD 17)

find_package(helib)

add_executable(bootstrapping main.cpp myutils.h)
target_link_libraries(bootstrapping helib ntl gmp gf2x)
target_include_directories(bootstrapping PUBLIC Helib/include)
cmake_minimum_required(VERSION 3.8)
project(aima)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp games.h)
add_executable(aima ${SOURCE_FILES})
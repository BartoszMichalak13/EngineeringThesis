cmake_minimum_required(VERSION 3.10)

# Project name
project(GraphGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Include the directory for header files
include_directories(include)

# Specify the executable target and its source files
add_executable(GraphGenerator
  src/main.cpp
  src/graph.cpp
  src/random.cpp
  src/utils.cpp
)
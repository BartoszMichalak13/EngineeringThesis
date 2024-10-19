cmake_minimum_required(VERSION 3.10)

# Project name
project(GraphGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Include the directory for header files
include_directories(include)

# Set optimization flag for Release mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -g") # as for 15.10.2024 almost x4 faster -O3 normally, without for debugining

# Specify the executable target and its source files
add_executable(GraphGenerator
  src/main.cpp
  src/graph.cpp
  src/random.cpp
  src/utils.cpp
)
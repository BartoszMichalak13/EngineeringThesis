cmake_minimum_required(VERSION 3.10)

# Project name
project(GraphGenerator)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Include the directory for header files
include_directories(include)


# Set optimization flag for Release mode
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -fsanitize=address -g -O3") # as for 15.10.2024 almost x4 faster -O3 normally, without for debugining
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -fsanitize=address -g") # as for 15.10.2024 almost x4 faster -O3 normally, without for debugining

# Specify the executable target and its source files
add_executable(GraphGenerator
  src/main.cpp
  src/graph.cpp
  src/random.cpp
  src/GraphStructures.cpp
  src/utils.cpp
  src/debugPrints.cpp
  src/helpers.cpp
  src/algorithmRunner.cpp
  src/Algorithms/TakahashiMatsuyama.cpp
  src/Algorithms/KouMarkowskyBerman.cpp
  src/Algorithms/DreyfusWagner.cpp
)
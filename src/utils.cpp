#include <iostream>
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"

// Function to parse input and update parameters
void parseInput(int argc, char *argv[], uint32_t &numberOfNodes, uint32_t &numberOfEdges, uint32_t &state)
{
  switch(argc) 
  {
    case 1: 
      std::cerr << "No arguments" << std::endl;
      state = 1;
      break;
    default:
      numberOfNodes = static_cast<uint32_t>(std::stoi(argv[1])); 
      numberOfEdges = static_cast<uint32_t>(std::stoi(argv[2])); 
      state = 0;
      break;
  }
}

// Function to output graph data or any relevant information
void writeOutput()
{
  // Placeholder for writing output to a file or console
  std::cout << "Output written successfully" << std::endl;
}

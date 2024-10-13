#include <iostream>
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"

int main(int argc, char *argv[]) 
{
  uint32_t numberOfNodes, numberOfEdges, state;

  //tmp solution
  float density = 0.1;

  parseInput(argc, argv, numberOfNodes, numberOfEdges, state);
  if(state) 
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }

  generateGraph(numberOfNodes, numberOfEdges, density);
  writeOutput();
  return 0;
}

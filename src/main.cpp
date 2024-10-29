#include <iostream>
// #include "graph.hpp"
#include "algorithmRunner.hpp"
// #include "random.hpp"
#include "utils.hpp"
#include <chrono>

int main(int argc, char *argv[]) 
{
  uint32_t numberOfNodes, numberOfEdges, state;

  //tmp solution
  float density = 0.1;
  bool printFlag = true;
  parseInput(argc, argv, numberOfNodes, numberOfEdges, state);
  if(state) 
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }
  auto start = std::chrono::high_resolution_clock::now();

  generateGraph(numberOfNodes, numberOfEdges, density, printFlag);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
  writeOutput();
  return 0;
}

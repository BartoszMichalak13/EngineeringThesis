#include <iostream>
#include "debugPrints.hpp"
#include "graph.hpp"
#include "algorithmRunner.hpp"
// #include "random.hpp"
#include "utils.hpp"
#include <chrono>
#include <memory>

int main(int argc, char *argv[]) 
{
  uint32_t numberOfNodes, numberOfEdges, state;
  std::vector<uint32_t> terminals(0);

  //tmp solution
  // float density = 0.1;
  // bool printFlag = true;
  bool printFlag = true;
  std::shared_ptr<Graph> graph = parseInput(argc, argv, numberOfNodes, numberOfEdges, state, printFlag, terminals);
  if(state) 
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }

  // std::shared_ptr<Graph> graph = generateGraph(numberOfNodes, numberOfEdges, printFlag);

  auto start = std::chrono::high_resolution_clock::now();

  runAlgorithms(graph, terminals);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
  writeOutput();

  return 0;
}

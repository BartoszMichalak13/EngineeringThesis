#include <iostream>
#include "debugPrints.hpp"
#include "graph.hpp"
#include "helpers.hpp"
#include "algorithmRunner.hpp"
// #include "random.hpp"
#include "utils.hpp"
#include <chrono>
#include <cmath>
#include <memory>

//TODO form of algorithmsToRun, density, loop to test, write output
//TODO REPAIR DREYFUS?? b04

int main(int argc, char *argv[]) 
{
  uint32_t numberOfNodes, numberOfEdges, state, algorithmsToRun, startingNumberOfNodes, targetNumberOfNodes;
  std::vector<uint32_t> terminals(0);

  std::string fileName;

  float density = 0.0;
  // bool printFlag = true;
  bool printFlag = false;
  // bool DreyfusPrintFlag = true;

  //TODO add file name to write to
  std::shared_ptr<Graph> graph = parseInput(
    argc,
    argv,
    fileName,
    numberOfNodes,
    numberOfEdges,
    density,
    algorithmsToRun,
    startingNumberOfNodes,
    targetNumberOfNodes,
    state,
    printFlag,
    terminals);

  if (state == 1000)
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }

  // std::shared_ptr<Graph> graph = generateGraph(numberOfNodes, numberOfEdges, printFlag);

  // auto start = std::chrono::high_resolution_clock::now();

  // check for 0 density?

  //TODO 1. make number of terminals more sensible
  //TODO 2. if not connected try creating new one, if idk 50 in row failed return eror CHECK IT
  //TODO 3. repair shortest path

  if (state == 2) {
    const uint32_t terminalOffset = 10;
    const uint32_t nodeOffset = 10;
    const uint32_t counterLimit = 50;
    uint32_t minNumberOfTerminals = std::ceil(startingNumberOfNodes / 50) + 2; // so there are minimum of 3
    for ( uint32_t numberOfTerminals = minNumberOfTerminals;
          numberOfTerminals < std::ceil(targetNumberOfNodes  * 0.9);
          numberOfTerminals += terminalOffset)
    {
      terminals = generateTerminals(startingNumberOfNodes, numberOfTerminals);
      for (; startingNumberOfNodes < targetNumberOfNodes; startingNumberOfNodes += nodeOffset) {

        //calculate number of Edges, gen graph
        uint32_t numberOfEdges = std::round(numberOfEdgesInClique(startingNumberOfNodes) * density);

        std::cout << "NumberOfNodes = " << startingNumberOfNodes << std::endl;
        std::cout << "numberOfEdges = " << numberOfEdges << std::endl;

        graph = generateGraph(startingNumberOfNodes, numberOfEdges, printFlag);
        std::pair<bool,bool> graphCheck = graph->isTree();
        uint32_t counter = 0;
        while (!graphCheck.second) {
          graph = generateGraph(startingNumberOfNodes, numberOfEdges, printFlag);
          graphCheck = graph->isTree();
          ++counter;
          if (counter >= counterLimit) {
            std::cerr << "Couldn't create connected graph" << std::endl;
            return 1;
          }
        }
        runAlgorithms(graph, terminals, algorithmsToRun, fileName);
      }
    }

  } else {
    std::cout << "HERE = " << algorithmsToRun << std::endl;

    runAlgorithms(graph, terminals, algorithmsToRun, fileName);
  }
  // auto stop = std::chrono::high_resolution_clock::now();
  // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  // std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;
  // writeOutput();

  return 0;
}

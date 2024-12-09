#include <iostream>
#include "debugPrints.hpp"
#include "graph.hpp"
#include "helpers.hpp"
#include "algorithmRunner.hpp"
#include "random.hpp"
#include "utils.hpp"
#include <chrono>
#include <cmath>
#include <memory>

int main(int argc, char *argv[])
{
  // Used for single instance generation
  uint32_t numberOfNodes, numberOfEdges;

  // Carries information about correctness and mode of operations
  uint32_t state;

  // Binary form of this number represents what algorithms should be run
  uint32_t algorithmsToRun;

  // Used as limits of node range in loop generation mode
  uint32_t LowerLimitOfNodes, UpperLimitOfNodes;

  // Terminal vertices which must be included in ST
  std::vector<uint32_t> terminals(0);

  // Destination to which we write our results
  std::string outputFileName;

  // Used in to generate graphs in mode 2 (loop generation)
  float density = 0.0;

  // Increases number of debug prints throughout the code
  bool printFlag = false;

  std::shared_ptr<Graph> graph = parseInput(
    argc,
    argv,
    outputFileName,
    numberOfNodes,
    numberOfEdges,
    density,
    algorithmsToRun,
    LowerLimitOfNodes,
    UpperLimitOfNodes,
    state,
    printFlag,
    terminals);

  // Checks on parsed input correctness
  if (state == 1000)
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }
  if (density <= 0 && state == 2) {
    std::cerr << "Non positive density " << density << std::endl;
    return state;
  }

  if (state == 2) // Loop generator
  {
    const uint32_t LowerLimitOfNodesCopy = LowerLimitOfNodes;

    //Modifies number of tests run
    const uint32_t numberOfTests = 15;
    for (uint32_t testNumber = 0; testNumber < numberOfTests; ++testNumber)
    {

      // Modifies growth rate of terminals in next tests
      const uint32_t terminalOffset = 1;

      const uint32_t maxNumberOfTerminals = 14;

      // Modifies growth rate of nodes in next tests
      const uint32_t nodeOffset = 10;

      // How many tris untill we decide that we have to few edges to make connected graph?
      const uint32_t counterLimit = 50;

      uint32_t minNumberOfTerminals = std::ceil(LowerLimitOfNodes / 50) + 2; // so there are minimum of 3

      // LowerLimitOfNodes = LowerLimitOfNodesCopy;
      for ( LowerLimitOfNodes = LowerLimitOfNodesCopy;
            LowerLimitOfNodes < UpperLimitOfNodes;
            LowerLimitOfNodes += nodeOffset)
      {
        // Calculate number of Edges to generate graph
        uint32_t numberOfEdges = std::round(numberOfEdgesInClique(LowerLimitOfNodes) * density);

        std::shared_ptr<Graph> g1 = generateGraph(LowerLimitOfNodes, numberOfEdges, printFlag);
        std::pair<bool,bool> graphCheck = g1->isTree();
        uint32_t counter = 0;

        // If graph instance was not connected try to re-generate it
        while (!graphCheck.second) {
          g1 = nullptr;
          g1 = generateGraph(LowerLimitOfNodes, numberOfEdges, printFlag);
          graphCheck = g1->isTree();
          ++counter;
          if (counter >= counterLimit)
          {
            std::cerr << "Couldn't create connected graph" << std::endl;
            return 1;
          }
        }

        for ( uint32_t numberOfTerminals = minNumberOfTerminals;
              numberOfTerminals < std::min(maxNumberOfTerminals,LowerLimitOfNodes);
              numberOfTerminals += terminalOffset)
        {

          std::vector<uint32_t> terminals1 = generateTerminals(LowerLimitOfNodes, numberOfTerminals);
          runAlgorithms(g1, terminals1, algorithmsToRun, outputFileName);
        }
        g1 = nullptr;
      }
    }
  }
  else // Either loaded graph or single generation mode
  {
    runAlgorithms(graph, terminals, algorithmsToRun, outputFileName);
  }

  // Clean singleton
  Random::destroyInstance();

  return 0;
}

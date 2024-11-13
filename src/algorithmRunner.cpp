#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "algorithmRunner.hpp"

//TODO split into 2 function
std::shared_ptr<Graph> generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool  printFlag)
{
  std::shared_ptr<Graph> graph(new Graph(numberOfNodes, numberOfEdges, printFlag));
  const uint32_t range = numberOfNodes - 1;

  for (uint32_t i = 0; i < numberOfEdges; ++i)
  {
    bool cannotMakeEdge = true;
    while(cannotMakeEdge) 
    {
      uint32_t node1Id = Random::getInstance()->generateRandomNumber(0, range);
      uint32_t node2Id = Random::getInstance()->generateRandomNumber(0, range);
      cannotMakeEdge = graph->checkIfEdgeExists(node1Id, node2Id);
      if (!cannotMakeEdge && node1Id != node2Id) //no multi edges, no self-loops
        graph->addEdgeWithRandomWeight(node1Id, node2Id);
    }
  }
  if (printFlag && numberOfNodes < 1000)
  {
    graph->printData();
  }
  return graph;
}

void runAlgorithms(std::shared_ptr<Graph> graph, std::vector<uint32_t> terminals)
{
  uint32_t numberOfNodes = graph->numberOfNodes;
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(graph, localCopyOfAdjacencyList);

  if (compareAdajcencyLists(graph->adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
    std::cout << "Comparing Adajcency Lists Failed" << std::endl;
    graph->printAdajcencyListFromGraph();
    std::cout << "Local" << std::endl;
    printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);
    return;
  }
  std::pair<bool,bool> graphCheck = graph->isTree();
  if (graphCheck.first) {
    std::cout << "Base Graph is acyclic" << std::endl;
  } else {
    std::cout << "Base Graph is NOT acyclic" << std::endl;
  }
  if (graphCheck.second) {
    std::cout << "Base Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  auto startMST = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> mst = graph->PrimMST(); //dummy for now

  auto stopMST = std::chrono::high_resolution_clock::now();
  auto durationMST = std::chrono::duration_cast<std::chrono::microseconds>(stopMST - startMST);

  graph->resetVisitedStatus();
  mst->resetVisitedStatus();

  std::cout << std::endl;
  std::pair<bool,bool> mstCheck = mst->isTree();
  if (mstCheck.first) {
    std::cout << "mst is acyclic" << std::endl;
    if (mstCheck.second) {
      std::cout << "mst is connected" << std::endl;
      std::cout << "mst is a Tree" << std::endl;
    } else {
      std::cout << "mst is NOT connected" << std::endl;
      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      return;
    }
  } else {
    std::cout << "mst is NOT acyclic" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  if (!terminals.size()) {
    uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
    terminals = graph->generateTerminals(numberOfTerminals);
  }

  if (graph->printFlag)
    printNodeVector(terminals); // print terminals

  auto startTMST = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> steinerTreeTakahashiMatsuyama = graph->TakahashiMatsuyama(terminals);

  auto stopTMST = std::chrono::high_resolution_clock::now();
  auto durationTMST = std::chrono::duration_cast<std::chrono::microseconds>(stopTMST - startTMST);

  graph->resetVisitedStatus();

  std::pair<bool,bool> steinerTreeTakahashiMatsuyamaCheck = steinerTreeTakahashiMatsuyama->isTree();
  if (steinerTreeTakahashiMatsuyamaCheck.first) {
    std::cout << "steinerTreeTakahashiMatsuyama is acyclic" << std::endl;
    if (steinerTreeTakahashiMatsuyamaCheck.second) {
      std::cout << "steinerTreeTakahashiMatsuyama is connected" << std::endl;
      std::cout << "steinerTreeTakahashiMatsuyama is a Tree" << std::endl;
    } else {
      std::cout << "steinerTreeTakahashiMatsuyama is NOT connected" << std::endl;
      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      return;
    }
  } else {
    std::cout << "steinerTreeTakahashiMatsuyama is NOT acyclic" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  if (compareAdajcencyLists(graph->adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
    std::cout << "Comparing Adajcency Lists Failed" << std::endl;
    graph->printAdajcencyListFromGraph();
    std::cout << "Local" << std::endl;
    printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  auto startKMBST = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> steinerTreeKouMarkowskyBerman(graph->KouMarkowskyBerman(terminals));

  auto stopKMBST = std::chrono::high_resolution_clock::now();
  auto durationKMBST = std::chrono::duration_cast<std::chrono::microseconds>(stopKMBST - startKMBST);
  graph->resetVisitedStatus();

  std::cout << std::endl;
  std::pair<bool,bool> steinerTreeKouMarkowskyBermanCheck = steinerTreeKouMarkowskyBerman->isTree();
  if (steinerTreeKouMarkowskyBermanCheck.first) {
    std::cout << "steinerTreeKouMarkowskyBerman is acyclic" << std::endl;
    if (steinerTreeKouMarkowskyBermanCheck.second) {
      std::cout << "steinerTreeKouMarkowskyBerman is connected" << std::endl;
      std::cout << "steinerTreeKouMarkowskyBerman is a Tree" << std::endl;
    } else {
      std::cout << "steinerTreeKouMarkowskyBerman is NOT connected" << std::endl;
      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      return;
    }
  } else {
    std::cout << "steinerTreeKouMarkowskyBerman is NOT acyclic" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << "COSTS:" << std::endl;
  // std::cout << "starting graph: " << graph->graphTotalCost() << std::endl;
  // std::cout << "MST: " << mst->graphTotalCost() << std::endl;
  // std::cout << "TakahashiMatsuyama steinerTree: " << steinerTreeTakahashiMatsuyama->graphTotalCost() << std::endl;
  // std::cout << "KouMarkowskyBerman steinerTree: " << steinerTreeKouMarkowskyBerman->graphTotalCost() << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << "Time:" << std::endl;
  // std::cout << "MST: " << durationMST.count() << " microseconds" << std::endl;
  // std::cout << "TakahashiMatsuyama steinerTree: " << durationTMST.count() << " microseconds" << std::endl;
  // std::cout << "KouMarkowskyBerman steinerTree: " << durationKMBST.count() << " microseconds" << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;





  uint32_t opt = graph->DreyfusWagner(terminals);

  graph->printAdajcencyListFromGraph();

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "COSTS:" << std::endl;
  std::cout << "starting graph: " << graph->graphTotalCost() << std::endl;
  std::cout << "MST: " << mst->graphTotalCost() << std::endl;
  std::cout << "TakahashiMatsuyama steinerTree: " << steinerTreeTakahashiMatsuyama->graphTotalCost() << std::endl;
  std::cout << "KouMarkowskyBerman steinerTree: " << steinerTreeKouMarkowskyBerman->graphTotalCost() << std::endl;
  std::cout << "DreyfusWagner(OPT) steinerTree: " << opt << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  // std::cout << "Time:" << std::endl;
  // std::cout << "MST: " << durationMST.count() << " microseconds" << std::endl;
  // std::cout << "TakahashiMatsuyama steinerTree: " << durationTMST.count() << " microseconds" << std::endl;
  // std::cout << "KouMarkowskyBerman steinerTree: " << durationKMBST.count() << " microseconds" << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;



  // graph->printAdajcencyListFromGraph();

  printNodeVector(terminals);

  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;
  return;
}

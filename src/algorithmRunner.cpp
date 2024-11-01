#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "algorithmRunner.hpp"


void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool  printFlag)
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
        graph->addEdge(node1Id, node2Id);
    }
  }
  if (printFlag && numberOfNodes < 1000)
  {
    graph->printData();
  }

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(graph, localCopyOfAdjacencyList);

  if (compareAdajcencyLists(graph->adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
    std::cout << "Comparing Adajcency Lists Failed" << std::endl;
    graph->printAdajcencyListFromGraph();
    std::cout << "Local" << std::endl;
    printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);
    return;
  }

  if(graph->isConnected()) {
    std::cout << "Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }



  std::shared_ptr<Graph> mst = graph->PrimMST(); //dummy for now

  graph->resetVisitedStatus();

  if(mst->isConnected()) { // TODO move code below to isConnected
    std::cout << "mst is connected" << std::endl;
  } else {
    std::cout << "mst is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
  std::vector<uint32_t> terminals = graph->generateTerminals(numberOfTerminals);

  printNodeVector(terminals);
  std::shared_ptr<Graph> steinerTreeTakahashiMatsuyama = graph->TakahashiMatsuyama(terminals);

  graph->resetVisitedStatus();

  if(steinerTreeTakahashiMatsuyama->isConnected()) { // TODO move code below to isConnected
    std::cout << "steinerTreeTakahashiMatsuyama is connected" << std::endl;
  } else {
    std::cout << "steinerTreeTakahashiMatsuyama is NOT connected" << std::endl;
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

  std::shared_ptr<Graph> steinerTreeKouMarkowskyBerman(graph->KouMarkowskyBerman(terminals));
  graph->resetVisitedStatus();

  if(steinerTreeKouMarkowskyBerman->isConnected()) { // TODO move code below to isConnected
    std::cout << "steinerTreeKouMarkowskyBerman is connected" << std::endl;
  } else {
    std::cout << "steinerTreeKouMarkowskyBerman is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }

  std::cout << "COSTS:" << std::endl;
  std::cout << "starting graph: " << graph->graphTotalCost() << std::endl;
  std::cout << "MST: " << mst->graphTotalCost() << std::endl;
  std::cout << "TakahashiMatsuyama steinerTree: " << steinerTreeTakahashiMatsuyama->graphTotalCost() << std::endl;
  std::cout << "KouMarkowskyBerman steinerTree: " << steinerTreeKouMarkowskyBerman->graphTotalCost() << std::endl;


  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;

}

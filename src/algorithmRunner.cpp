#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "algorithmRunner.hpp"


void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool  printFlag)
{
  std::shared_ptr<Graph> graph(new Graph(numberOfNodes, numberOfEdges, printFlag));
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraph(graph, localCopyOfAdjacencyList);
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;


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

  if(graph->isConnected()) {
    std::cout << "Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;

  std::shared_ptr<Graph> mst = graph->PrimMST(); //dummy for now
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;

  graph->resetVisitedStatus();

  if(mst->isConnected()) { // TODO move code below to isConnected
    std::cout << "mst is connected" << std::endl;
  } else {
    std::cout << "mst is NOT connected" << std::endl;
    delete[] localCopyOfAdjacencyList;
    localCopyOfAdjacencyList = nullptr;
    return;
  }
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;

  uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
  std::vector<uint32_t> terminals = graph->generateTerminals(numberOfTerminals);

  // graph->printAdajcencyListFromGraph();
  std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;
  printNodeVector(terminals);
  std::shared_ptr<Graph> steinerTreeTakahashiMatsuyama = graph->TakahashiMatsuyama(terminals);
    std::cout << "HELLO" << std::endl;
  if (graph == nullptr) {
    std::cout << "nullptr" << std::endl;
  } else {
    std::cout << "ok" << std::endl;
    if (graph->numberOfNodes == 0) {
      std::cout << "NULL" << std::endl;
    } else {
      std::cout << "graph->numberOfNodes " << graph->numberOfNodes << std::endl;
      std::cout << "steinerTreeTakahashiMatsuyama->numberOfNodes " << steinerTreeTakahashiMatsuyama->numberOfNodes << std::endl;
    }
  }

  graph->printAdajcencyListFromGraph();

  graph->resetVisitedStatus();
  std::cout << "HELLO" << std::endl;

  if(steinerTreeTakahashiMatsuyama->isConnected()) { // TODO move code below to isConnected
    std::cout << "steinerTreeTakahashiMatsuyama is connected" << std::endl;
  } else {
    std::cout << "steinerTreeTakahashiMatsuyama is NOT connected" << std::endl;
    // return;
  }

  std::cout << "HELLO" << std::endl;

  if (compareAdajcencyLists(graph->adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
    std::cout << "Comparing Adajcency Lists Failed" << std::endl;
    // return dummySharedPointerGraph();
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
  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;

}

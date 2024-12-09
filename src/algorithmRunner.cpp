#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "algorithmRunner.hpp"

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

void runAlgorithms(std::shared_ptr<Graph> graph, std::vector<uint32_t> terminals, uint32_t algorithmsToRun, const std::string& fileName)
{
  uint32_t numberOfNodes = graph->numberOfNodes;
  bool pF = graph->printFlag;

  std::pair<bool,bool> graphCheck = graph->isTree();
  if (graphCheck.first) {
    if (pF) std::cout << "Base Graph is acyclic" << std::endl;
  } else {
    if (pF) std::cout << "Base Graph is NOT acyclic" << std::endl;
  }
  if (graphCheck.second) {
    if (pF) std::cout << "Base Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    return;
  }
  /*******************************************************************/
  /*                             MST                                 */
  /*******************************************************************/
  auto startMST = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> mst = graph->PrimMST(); //dummy for now

  auto stopMST = std::chrono::high_resolution_clock::now();
  auto durationMST = std::chrono::duration_cast<std::chrono::microseconds>(stopMST - startMST);

  graph->resetVisitedStatus();
  mst->resetVisitedStatus();

  std::cout << std::endl;
  std::pair<bool,bool> mstCheck = mst->isTree();
  if (mstCheck.first) {
    if (pF) std::cout << "mst is acyclic" << std::endl;
    if (mstCheck.second) {
      if (pF) std::cout << "mst is connected" << std::endl;
      if (pF) std::cout << "mst is a Tree" << std::endl;
    } else {
      std::cout << "mst is NOT connected" << std::endl;
      return;
    }
  } else {
    std::cout << "mst is NOT acyclic" << std::endl;
    return;
  }

  /*            GENERATE TERMINALS IF NOT GIVEN (MODE 1)             */
  if (!terminals.size()) {
    uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
    terminals = graph->generateTerminals(numberOfTerminals);
  }

  if (pF)
    printUintVector(terminals); // print terminals


  /*******************************************************************/
  /*                      Takahashi Matsuyama                        */
  /*******************************************************************/
  uint32_t TMSTcost = 0;
  std::chrono::microseconds durationTMST = std::chrono::microseconds(0);
  /*
    duration of not necessary
    duration of init
    duration of preparation for next iteration
    duration of main loop
  */
  std::vector<std::chrono::microseconds> TMtimes = std::vector<std::chrono::microseconds>(10, std::chrono::microseconds(0));

  if (isNthBitSet(algorithmsToRun, 0)) {

    auto startTMST = std::chrono::high_resolution_clock::now();

    std::shared_ptr<Graph> steinerTreeTakahashiMatsuyama = graph->TakahashiMatsuyama(terminals, TMtimes);

    auto stopTMST = std::chrono::high_resolution_clock::now();
    durationTMST = std::chrono::duration_cast<std::chrono::microseconds>(stopTMST - startTMST);
    TMSTcost = steinerTreeTakahashiMatsuyama->graphTotalCost();
    graph->resetVisitedStatus();

    std::pair<bool,bool> steinerTreeTakahashiMatsuyamaCheck = steinerTreeTakahashiMatsuyama->isTree();
    if (steinerTreeTakahashiMatsuyamaCheck.first) {
      if(pF) std::cout << "steinerTreeTakahashiMatsuyama is acyclic" << std::endl;
      if (steinerTreeTakahashiMatsuyamaCheck.second) {
        if(pF) std::cout << "steinerTreeTakahashiMatsuyama is connected" << std::endl;
        if(pF) std::cout << "steinerTreeTakahashiMatsuyama is a Tree" << std::endl;
      } else {
        std::cout << "steinerTreeTakahashiMatsuyama is NOT connected" << std::endl;
        return;
      }
    } else {
      std::cout << "steinerTreeTakahashiMatsuyama is NOT acyclic" << std::endl;
      return;
    }
  }

  /*******************************************************************/
  /*                     KOU MARKOWSKI BERMAN                        */
  /*******************************************************************/
  uint32_t KMBSTcost = 0;
  std::chrono::microseconds durationKMBST = std::chrono::microseconds(0);
  /*
    duration of not necessary
    duration of step1
    duration of step2
    duration of step3
    duration of step4
    duration of step5
  */
  std::vector<std::chrono::microseconds> KMBtimes = std::vector<std::chrono::microseconds>(10, std::chrono::microseconds(0));
  if (isNthBitSet(algorithmsToRun, 1)) {
    auto startKMBST = std::chrono::high_resolution_clock::now();

    std::shared_ptr<Graph> steinerTreeKouMarkowskyBerman(graph->KouMarkowskyBerman(terminals, KMBtimes));

    auto stopKMBST = std::chrono::high_resolution_clock::now();
    durationKMBST = std::chrono::duration_cast<std::chrono::microseconds>(stopKMBST - startKMBST);
    KMBSTcost = steinerTreeKouMarkowskyBerman->graphTotalCost();
    graph->resetVisitedStatus();

    std::cout << std::endl;
    std::pair<bool,bool> steinerTreeKouMarkowskyBermanCheck = steinerTreeKouMarkowskyBerman->isTree();
    if (steinerTreeKouMarkowskyBermanCheck.first) {
      if(pF) std::cout << "steinerTreeKouMarkowskyBerman is acyclic" << std::endl;
      if (steinerTreeKouMarkowskyBermanCheck.second) {
        if(pF) std::cout << "steinerTreeKouMarkowskyBerman is connected" << std::endl;
        if(pF) std::cout << "steinerTreeKouMarkowskyBerman is a Tree" << std::endl;
      } else {
        std::cout << "steinerTreeKouMarkowskyBerman is NOT connected" << std::endl;
        return;
      }
    } else {
      std::cout << "steinerTreeKouMarkowskyBerman is NOT acyclic" << std::endl;
      return;
    }
  }

  /*******************************************************************/
  /*                         DREYFUS WAGNER                          */
  /*******************************************************************/
  uint32_t DWSTcost = 0;
  std::chrono::microseconds durationDreyfus = std::chrono::microseconds(0);
  if (isNthBitSet(algorithmsToRun, 2)) {
    auto startDreyfus = std::chrono::high_resolution_clock::now();

    DWSTcost = graph->DreyfusWagner(terminals);

    auto stopDreyfus = std::chrono::high_resolution_clock::now();
    durationDreyfus = std::chrono::duration_cast<std::chrono::microseconds>(stopDreyfus - startDreyfus);
  }


  if(pF) {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "COSTS:" << std::endl;
    std::cout << "starting graph: " << graph->graphTotalCost() << std::endl;
    std::cout << "MST: " << mst->graphTotalCost() << std::endl;

    if (isNthBitSet(algorithmsToRun, 0))
      std::cout << "TakahashiMatsuyama steinerTree: " << TMSTcost << std::endl;

    if (isNthBitSet(algorithmsToRun, 1))
      std::cout << "KouMarkowskyBerman steinerTree: " << KMBSTcost << std::endl;

    if (isNthBitSet(algorithmsToRun, 2))
      std::cout << "DreyfusWagner(OPT) steinerTree: " << DWSTcost << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Time:" << std::endl;
    std::cout << "MST: " << durationMST.count() << " microseconds" << std::endl;

    if (isNthBitSet(algorithmsToRun, 0))
      std::cout << "TakahashiMatsuyama steinerTree: " << durationTMST.count() << " microseconds" << std::endl;

    if (isNthBitSet(algorithmsToRun, 1))
      std::cout << "KouMarkowskyBerman steinerTree: " << durationKMBST.count() << " microseconds" << std::endl;

    if (isNthBitSet(algorithmsToRun, 2))
      std::cout << "DreyfusWagner(OPT) steinerTree: " << durationDreyfus.count() << " microseconds" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }

  if (pF)
    printUintVector(terminals);

  if (fileName != "") {
    writeOutput(
      fileName,
      terminals.size(),
      graph->numberOfNodes,
      graph->numberOfEdges,
      DWSTcost,
      durationDreyfus.count(),
      TMSTcost,
      durationTMST.count(),
      TMtimes,
      KMBSTcost,
      durationKMBST.count(),
      KMBtimes
    );
  }

  return;
}

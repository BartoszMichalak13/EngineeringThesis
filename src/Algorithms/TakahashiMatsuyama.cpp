#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"

/*
  As Takahashi and Matsuyama state in their article, this algorithm has 2 key steps.
  * 1.*
      Add single terminal to the Tree T_0
  * 2.*
      for i in range (0,n-1)
        search for the nearest terminal to the tree
        add it with shortest path leading to it to T_i creating T_i+1
*/
std::shared_ptr<Graph> Graph::TakahashiMatsuyama(
    std::vector<uint32_t> terminals,
    std::vector<std::chrono::microseconds> &timeMeasurements)
{
  /*                                              Time Measurement Variables                                                     */
  std::chrono::time_point<std::chrono::high_resolution_clock> startNotNecessary{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopNotNecessary{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startPrepareForNextIteration{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopPrepareForNextIteration{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startInit{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopInit{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startMainLoop{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopMainLoop{std::chrono::high_resolution_clock::duration{0}};

  std::chrono::microseconds durationNotNecessary{0};
  std::chrono::microseconds durationPrepareForNextIteration{0};
  std::chrono::microseconds durationInit{0};
  std::chrono::microseconds durationMainLoop{0};

  /***                 Time for all copies which are not necessary by algorithm itself to work properly                        ***/
  startNotNecessary = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> self = shared_from_this();

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  resetVisitedStatus();
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(self, localCopyOfAdjacencyList);

  stopNotNecessary = std::chrono::high_resolution_clock::now();
  durationNotNecessary += std::chrono::duration_cast<std::chrono::microseconds>(stopNotNecessary - startNotNecessary);

  std::vector<PseudoEdge> tmpTreeEdges;

  bool foundTerminal = false;
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;

  /***                 Time for all copies which are not necessary by algorithm itself to work properly                        ***/
  startInit = std::chrono::high_resolution_clock::now();

  /*                                  Insert all neighbours of first terminal to queue                                           */
  std::shared_ptr<Edge> selfLoopInitEdge = constructSelfLoopInitEdge(vertices[terminals.at(0)]);
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[terminals.at(0)].size(); ++i)
  {
    updatePred(localCopyOfAdjacencyList[terminals.at(0)].at(i), selfLoopInitEdge, localCopyOfAdjacencyList);
    toVisit.push(localCopyOfAdjacencyList[terminals.at(0)].at(i));
  }

  vertices[terminals.at(0)]->visited = true;
  terminals.erase(terminals.begin());

  stopInit = std::chrono::high_resolution_clock::now();
  durationInit += std::chrono::duration_cast<std::chrono::microseconds>(stopInit - startInit);

  do
  {
    /*                          Empty queue and repopulate it with edges from original graph and 0-edges                         */
    /* Reason: As edge weights are changed during process of finding shortest path, we then make this shortest path part of a ST */
    /* That means - during the process of finding next shortest path, we measure distance from the tree so obviously distances   */
    /* from freshly added edges are lower than in the queue (since they were distances from previous, smaller tree). So we need  */
    /* to reset them.                                                                                                            */
    if (foundTerminal)
    {
      startPrepareForNextIteration = std::chrono::high_resolution_clock::now();

      // reset priority queue
      toVisit = std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers>();

      /*                             Add neighbours of all unique nodes in tree to the queue                                     */
      std::vector<uint32_t> uniqueNodes = tmpPseudoEdgeReturnUniqueNumbers(tmpTreeEdges);//TODO test it
      for (uint32_t i = 0; i < uniqueNodes.size(); ++i)
      {
        int32_t idx = findInArray(uniqueNodes.at(i), vertices, numberOfNodes);
        if (idx == -1) {
          std::cerr << "Error: Node " << uniqueNodes.at(i) << " not found in queue reset and repopulation"  << std::endl;
          return dummySharedPointerGraph();
        }

        std::shared_ptr<Edge> e = findZeroEdgeInAdjacentTo(idx, localCopyOfAdjacencyList);
        if (e == nullptr) {
          std::cerr << "Error: e doesnt exit aka no neighbour with 0 weight for node "<< uniqueNodes.at(i) << std::endl;
          return dummySharedPointerGraph();
        }

        searchNeighboursV2(toVisit, localCopyOfAdjacencyList, uniqueNodes.at(i), e);
      }
      foundTerminal = false;

      stopPrepareForNextIteration = std::chrono::high_resolution_clock::now();
      durationPrepareForNextIteration += std::chrono::duration_cast<std::chrono::microseconds>(stopPrepareForNextIteration - startPrepareForNextIteration);
    }

    /*        MAIN LOOP - here we add nodes to queue and check if they are terminals if so we augment the current tree.          */
    while(!foundTerminal)
    {
      startMainLoop = std::chrono::high_resolution_clock::now();
      // remoeve all edges that lead to already visited nodes
      if (!toVisit.empty())
        while (toVisit.top() == nullptr || toVisit.top()->end->visited) //TODO is it always the end?
          toVisit.pop();

      if (toVisit.empty()) {
        std::cerr << "Error: End in loop TakahashiMatsuyama; dummy graph" << std::endl;
        return dummySharedPointerGraph();
      }

      std::shared_ptr<Edge> e = toVisit.top();
      toVisit.pop();
      e->end->visited = true;
      uint32_t nextNodeIndex = e->end->id;

      // check if we have found terminal
      int32_t idx = findInUintVector(nextNodeIndex, terminals);
      if (idx > -1) {
        foundTerminal = true;
        // remove found terminal from the list
        terminals.erase(terminals.begin() + idx);
        std::shared_ptr<Edge> oldPred = e->pred;
        embedEdgeIntoTree(e, localCopyOfAdjacencyList, adjacencyList, tmpTreeEdges);

        // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
        while (oldPred != nullptr && oldPred->start->id != e->start->id) {
          e = oldPred;
          oldPred = e->pred;
          embedEdgeIntoTree(e, localCopyOfAdjacencyList, adjacencyList, tmpTreeEdges);
        } // untill we reach beginning of the path
        resetVisitedStatusAndWeightInCopyOfAdjacencyList(localCopyOfAdjacencyList, self); // reset visitied status
      }
      else // aka terminal not found
      {
        searchNeighboursV2(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
      }
      stopMainLoop = std::chrono::high_resolution_clock::now();
      durationMainLoop += std::chrono::duration_cast<std::chrono::microseconds>(stopMainLoop - startMainLoop);
    }
  } while(!terminals.empty());

  uint64_t totalWeight = 0;
  for (uint32_t i = 0; i < tmpTreeEdges.size(); ++i)
    totalWeight += tmpTreeEdges.at(i).weight;

  /***                   Time for all the prints and creating tree structure (only value would suffice)                        ***/
  startNotNecessary = std::chrono::high_resolution_clock::now();

  if (printFlag) {
    std::cout << "TakahashiMatsuyama totalWeight = " << totalWeight << std::endl;
  }
  for (uint32_t i = 0; i < originalTerminals.size(); ++i) {
    bool found = false;
    for (uint32_t j = 0; j < tmpTreeEdges.size(); ++j) {
      if (originalTerminals.at(i) == tmpTreeEdges.at(j).start || originalTerminals.at(i) == tmpTreeEdges.at(j).end) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << std::endl;
      std::cerr << "Error: missing originalTerminals: " << originalTerminals.at(i) << std::endl;
    }
  }
  if (printFlag) std::cout << std::endl;

  std::vector<uint32_t> nodes = tmpPseudoEdgeReturnUniqueNumbers(tmpTreeEdges);

  if (printFlag && numberOfNodes < 1000) {
    std::cout << "Edges in  tmpTreeEdges: " << std::endl;
    for (uint32_t i = 0; i < tmpTreeEdges.size(); ++i)
      std::cout << tmpTreeEdges.at(i).start << "->" <<  tmpTreeEdges.at(i).end << "; ";
    std::cout << std::endl;
    std::cout << "Original terminals: ";
    for (uint32_t i = 0; i < originalTerminals.size(); ++i)
      std::cout << originalTerminals.at(i) << ", ";
    std::cout << std::endl;
    std::cout << "Nodes in Tree: ";
    for (uint32_t i = 0; i < nodes.size(); ++i)
      std::cout << nodes.at(i) << ", ";
    std::cout << std::endl;
  }

  std::vector<std::shared_ptr<Edge>> treeEdges;
  for (uint32_t i = 0; i < tmpTreeEdges.size(); ++i) {
    std::shared_ptr<Edge> e = findEdge(tmpTreeEdges.at(i).start, tmpTreeEdges.at(i).end, adjacencyList);
    if (e == nullptr) {
      std::cerr << "Error: missing tmpTreeEdges: " << tmpTreeEdges.at(i).start << " - " << tmpTreeEdges.at(i).end << std::endl;
      return dummySharedPointerGraph();
    } else{
      treeEdges.push_back(e);
    }
  }

  /*                                            Steiner Tree creation                                                            */
  std::shared_ptr<Graph> steinerTree(new Graph(nodes, treeEdges, nodes.size(), treeEdges.size(), printFlag));

  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    for (uint32_t j = 0; j < localCopyOfAdjacencyList[i].size(); ++j) {
      localCopyOfAdjacencyList[i].at(j)->start = nullptr;
      localCopyOfAdjacencyList[i].at(j)->end = nullptr;
      localCopyOfAdjacencyList[i].at(j)->pred = nullptr;
      localCopyOfAdjacencyList[i].at(j)->succ = nullptr;
    }
    localCopyOfAdjacencyList[i].clear();
  }

  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;
  resetVisitedStatus();

  stopNotNecessary = std::chrono::high_resolution_clock::now();
  durationNotNecessary += std::chrono::duration_cast<std::chrono::microseconds>(stopNotNecessary - startNotNecessary);

  timeMeasurements.push_back(durationNotNecessary);
  timeMeasurements.push_back(durationInit);
  timeMeasurements.push_back(durationPrepareForNextIteration);
  timeMeasurements.push_back(durationMainLoop);
  return steinerTree;
}
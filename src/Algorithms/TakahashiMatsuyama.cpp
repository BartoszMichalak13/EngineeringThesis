#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"

//TODO repair, sometimes not acyclic
/*
Edges in  tmpTreeEdges:
5->9; 12->5; 0->12; 7->0; 5->0; 9->8;
*/


/*difference between prim -
1. we update weights in edges, based on distance to beg. node
2. in this variation of TakahashiMatsuyama we search just untill we meet node we search for
3. Terminals should be narrowed down after each found one
KIEDY ZNAJDE WAGI W NOWYM GRAFIE NA 0000000000000000000

usuwam wiekszosc pointerow, gdyz tworze nowe drzewo, a nie chce uszkadzac starego

modify edges, give pointer to pred

@param terminals should be of positive size
*/
std::shared_ptr<Graph> Graph::TakahashiMatsuyama(std::vector<uint32_t> terminals) {// do it with priority Queue, maybe my own immplementation?
  std::shared_ptr<Graph> self = shared_from_this();

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  //init
  resetVisitedStatus();
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(self, localCopyOfAdjacencyList);


  //TODO change it to normal Edges
  std::vector<PseudoEdge> tmpTreeEdges;

  bool foundTerminal = false;
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;

  //TODO make it safe
  std::shared_ptr<Edge> selfLoopInitEdge = constructSelfLoopInitEdge(vertices[terminals.at(0)]);

  for (uint32_t i = 0; i < localCopyOfAdjacencyList[terminals.at(0)].size(); ++i) {
    updatePred(localCopyOfAdjacencyList[terminals.at(0)].at(i), selfLoopInitEdge, localCopyOfAdjacencyList);
    toVisit.push(localCopyOfAdjacencyList[terminals.at(0)].at(i));
  }

  vertices[terminals.at(0)]->visited = true;
  terminals.erase(terminals.begin());

  do
  {
    if (foundTerminal)
    {
      while(!toVisit.empty())
        toVisit.pop();

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
    }

    //main loop
    while(!foundTerminal)
    { // O(n^2)
      //remoeve all edges that lead to already visited nodes
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
        //remove found terminal from the list
        terminals.erase(terminals.begin() + idx);
        std::shared_ptr<Edge> oldPred = e->pred;
        embedEdgeIntoTree(e, localCopyOfAdjacencyList, adjacencyList, tmpTreeEdges);

        // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
        while (oldPred != nullptr && oldPred->start->id != e->start->id) {
          e = oldPred;
          oldPred = e->pred;
          embedEdgeIntoTree(e, localCopyOfAdjacencyList, adjacencyList, tmpTreeEdges);
        } // untill we reach beginning of the path
        resetVisitedStatusInCopyOfAdjacencyList(localCopyOfAdjacencyList, self); //it should reset visitied status
      } else { // aka terminal not found
        searchNeighboursV2(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
      }
    }
  } while(!terminals.empty()); // O(k)

  uint64_t totalWeight = 0;
  for (uint32_t i = 0; i < tmpTreeEdges.size(); ++i)
    totalWeight += tmpTreeEdges.at(i).weight;

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
  std::cout << std::endl;
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
  // std::cout << "numberOfNodes " << numberOfNodes << "; steinerTree->numberOfNodes " << steinerTree->numberOfNodes << std::endl;
  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;
  resetVisitedStatus();
  return steinerTree;
}
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <algorithm>

//TODO na spokojnie sprawdz czy all pointery dzialaja tj, czy sa dobrze przepisywane
/*
1. Construct complete graph G_1 where V = terminals (Steiner Points), E = shortest paths between terminals
2. Find mst of G_1
3. Construct graph G_s by replacing each edge from mst with corresponding shortest path
4. Find mst of G_s
5. Remove all branches that lead to leaf which is not Steiner Point (terminal)

In all steps above if there are multiple answers to choose from (e.g. multiple msts) choos aribitrary one

TODO Possible imporvements: what if we create not a clique, but a "normal" graph, we use same method, but we do not contract edges
We may check basic algorithm and my changed version
Most likely (almost obvoiusly) it will only change execution time, not resulting weight of Steiner Tree
*/
std::shared_ptr<Graph> Graph::KouMarkowskyBerman(std::vector<uint32_t> terminals) {

  const uint32_t numberOfSteinerPoints = terminals.size();
  const uint32_t numberOfCliqueEdges = numberOfEdgesInClique(numberOfSteinerPoints);

  printNodeVector(terminals);

  std::vector<std::shared_ptr<Edge>> ShortestPaths;
  std::shared_ptr<std::vector<std::shared_ptr<Edge>>>* tmpShortestPaths = new std::shared_ptr<std::vector<std::shared_ptr<Edge>>>[numberOfCliqueEdges];

  uint32_t currentShortestPathIdx = 0;
  for (uint32_t i = 0; i < numberOfSteinerPoints; ++i) {
    for (uint32_t j = i + 1; j < numberOfSteinerPoints; ++j) {
      std::cout << "ShortestPath(" << terminals.at(i) << ", " << terminals.at(j) << ");" << std::endl;
      tmpShortestPaths[currentShortestPathIdx] = ShortestPath(terminals.at(i), terminals.at(j));
      ++currentShortestPathIdx;
    }
  }

  for (uint32_t i = 0; i < numberOfCliqueEdges; ++i) {
    uint32_t shortestPathWeight = 0;
    for (uint32_t j = 0; j < tmpShortestPaths[i]->size(); ++j) {
      shortestPathWeight += tmpShortestPaths[i]->at(j)->weight;
    }
    std::shared_ptr<Node> start = tmpShortestPaths[i]->at(0)->start;
    std::shared_ptr<Node> end = tmpShortestPaths[i]->at(tmpShortestPaths[i]->size() - 1)->end;
    std::cout << "start " << tmpShortestPaths[i]->at(0)->start->id << "; end " << tmpShortestPaths[i]->at(tmpShortestPaths[i]->size() - 1)->end->id << std::endl;
    std::cout << "alt start " << tmpShortestPaths[i]->at(0)->end->id << "; alt end " << tmpShortestPaths[i]->at(tmpShortestPaths[i]->size() - 1)->start->id << std::endl;
    ShortestPaths.push_back(std::shared_ptr<Edge>(new Edge(start, shortestPathWeight, end)));
  }
// 30, 24, 14, 32, 22, 26, 10, 34, 0, 21,
  uint32_t counter = 0;
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  printEdgeVector(ShortestPaths);
  printNodeVector(terminals);
  std::shared_ptr<Graph> g1(new Graph(terminals, ShortestPaths, numberOfSteinerPoints, numberOfCliqueEdges, printFlag));
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  std::shared_ptr<Graph> t1 = g1->PrimMST();
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  std::vector<uint32_t> treeNodes;
  std::vector<std::shared_ptr<Edge>> treeEdges;
  uint32_t numberOfTreeNodes = 0;
  uint32_t numberOfTreeEdges = 0;
  //przeiteruj cale adj.list t1 i zamien na tmpShortestPaths
  //Trzeba wyekstrachowac wierzcholki i krawędzie
  // BEZ POWTÓRZEN KRAW I WIERZCH
  for(uint32_t i = 0; i < numberOfSteinerPoints; ++i) {
    for (uint32_t j = 0; j < t1->adjacencyList[i].size(); ++j) {
      std::shared_ptr<Node> start = t1->adjacencyList[i].at(j)->start;
      std::shared_ptr<Node> end = t1->adjacencyList[i].at(j)->end;
      int32_t idx = findInEdgeVector(start, end, ShortestPaths);
      if (idx > -1) {
        for (uint32_t k = 0; k < tmpShortestPaths[idx]->size(); ++k) {
          //TODO make part below more transparent possibly add addEdgeIfNotAlreadyIn
          //we can do that bc ShortestPaths.at(i) corresponds to tmpShortestPaths[i]
          std::shared_ptr<Edge> e = tmpShortestPaths[idx]->at(k);
          int32_t repetitionCheck = findInEdgeVector(e->start, e->end, treeEdges);
          if (repetitionCheck == -1) {
            treeEdges.push_back(e);
            ++numberOfTreeEdges;

            addNodeIfNotAlreadyIn(e->start->id, treeNodes, numberOfTreeNodes);
            addNodeIfNotAlreadyIn(e->end->id, treeNodes, numberOfTreeNodes);
          }
        }
      }
    }
  }
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  std::shared_ptr<Graph> gs(new Graph(treeNodes, treeEdges, numberOfTreeNodes, numberOfTreeEdges, printFlag));
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  std::shared_ptr<Graph> ts = gs->PrimMST();
  std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  uint32_t numberOfSteinerTreeNodes     = numberOfTreeNodes;
  uint32_t tmpNumberOfSteinerTreeNodes  = numberOfTreeNodes;

  uint32_t numberOfSteinerTreeEdges     = numberOfTreeEdges;
  uint32_t tmpNumberOfSteinerTreeEdges  = numberOfTreeEdges;

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfTreeNodes];
  copyAdjacencyListFromGraph(ts, localCopyOfAdjacencyList);

  bool foundDanglindBnrach = false;
  bool continueSearch = true;
  //TODO delete leaves
  while (continueSearch) {
    for (uint32_t i = 0; i < numberOfSteinerTreeNodes; ++i) {
      //if non terminal leaf

      //TODO check if it's safe (e.i. if start == localCopyOfAdjacencyList[i])
      if (localCopyOfAdjacencyList[i].size() == 1 && findInUintVector(localCopyOfAdjacencyList[i].at(0)->start->id, terminals) == -1) {
        foundDanglindBnrach = true;

        treeNodes.erase(std::remove(treeNodes.begin(), treeNodes.end(), localCopyOfAdjacencyList[i].at(0)->start->id), treeNodes.end());
        --tmpNumberOfSteinerTreeNodes;

        std::shared_ptr<Edge> e = findInEdgeVectorAndReturnValue(localCopyOfAdjacencyList[i].at(0)->start, localCopyOfAdjacencyList[i].at(0)->end, treeEdges);
        if (e != nullptr){
          treeEdges.erase(std::remove(treeEdges.begin(), treeEdges.end(), e), treeEdges.end());
          --tmpNumberOfSteinerTreeEdges;
          removeEdgeFromAdjacencyList(e, localCopyOfAdjacencyList);
        } else {
          std::cerr << "Edge not found in delete leaes part" << std::endl;
          return dummySharedPointerGraph();
        }
      }
    }
    if (!foundDanglindBnrach) {
      continueSearch = false;
      // break;
    } else {
      foundDanglindBnrach = false;
      sortAdjacencyList(numberOfSteinerTreeNodes, localCopyOfAdjacencyList);
      numberOfSteinerTreeNodes = tmpNumberOfSteinerTreeNodes;
      numberOfSteinerTreeEdges = tmpNumberOfSteinerTreeEdges;
    }
  }
    std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  resetVisitedStatus();
    std::cout << "steinerTreeKouMarkowskyBerman "<< counter++ << std::endl;

  std::shared_ptr<Graph> th(new Graph(treeNodes, treeEdges, numberOfSteinerTreeNodes, numberOfSteinerTreeEdges, printFlag));
  return th;
}
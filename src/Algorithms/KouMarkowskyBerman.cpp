#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <algorithm>

//TODO ASDFASDFASFASFASFA
// 7->10; 7->11; 10->7; 10->11; 11->7; 11->5; 11->10; 11->1; 11->3; 5->11; 1->11; 3->11; 3->12; 12->3;
// steinerTreeKouMarkowskyBerman is NOT acyclic

//TODO na spokojnie sprawdz czy all pointery dzialaja tj, czy sa dobrze przepisywane
//TODO opis
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
std::shared_ptr<Graph> Graph::KouMarkowskyBerman(
    std::vector<uint32_t> terminals,
    std::vector<std::chrono::microseconds> &timeMeasurements)
{
  /*                                              Time Measurement Variables                                                     */
  std::chrono::time_point<std::chrono::high_resolution_clock> startNotNecessary{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopNotNecessary{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startStep1{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopStep1{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startStep2{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopStep2{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startStep3{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopStep3{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startStep4{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopStep4{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> startStep5{std::chrono::high_resolution_clock::duration{0}};
  std::chrono::time_point<std::chrono::high_resolution_clock> stopStep5{std::chrono::high_resolution_clock::duration{0}};

  std::chrono::microseconds durationNotNecessary{0};
  std::chrono::microseconds durationStep1{0};
  std::chrono::microseconds durationStep2{0};
  std::chrono::microseconds durationStep3{0};
  std::chrono::microseconds durationStep4{0};
  std::chrono::microseconds durationStep5{0};

  startNotNecessary = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> self = shared_from_this();

  stopNotNecessary = std::chrono::high_resolution_clock::now();
  durationNotNecessary += std::chrono::duration_cast<std::chrono::microseconds>(stopNotNecessary - startNotNecessary);

  /******************************************************************************************************************************/
  /*                         Step 1 - Create CLique on terminals from shortest paths between them                               */
  /******************************************************************************************************************************/
  startStep1 = std::chrono::high_resolution_clock::now();

  const uint32_t numberOfSteinerPoints = terminals.size();
  const uint32_t numberOfCliqueEdges = numberOfEdgesInClique(numberOfSteinerPoints);
  
  std::vector<std::shared_ptr<Edge>> ShortestPaths;
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> tmpShortestPaths = AllPairsShortestPath(terminals);// new std::shared_ptr<std::vector<std::shared_ptr<Edge>>>[numberOfCliqueEdges];

  for (uint32_t i = 0; i < numberOfCliqueEdges; ++i) {
    uint32_t shortestPathWeight = 0;
    for (uint32_t j = 0; j < tmpShortestPaths.at(i)->size(); ++j) {
      shortestPathWeight += tmpShortestPaths.at(i)->at(j)->weight;
    }
    // TODO this comment:
    // tmpShortestPaths is "kinda" inverted
    std::shared_ptr<Node> start = tmpShortestPaths.at(i)->at(0)->end;
    std::shared_ptr<Node> end = tmpShortestPaths.at(i)->at(tmpShortestPaths.at(i)->size() - 1)->start;
    ShortestPaths.push_back(std::shared_ptr<Edge>(new Edge(start, shortestPathWeight, end)));
  }
  std::shared_ptr<Graph> g1(new Graph(terminals, ShortestPaths, numberOfSteinerPoints, numberOfCliqueEdges, printFlag));

  stopStep1 = std::chrono::high_resolution_clock::now();
  durationStep1 += std::chrono::duration_cast<std::chrono::microseconds>(stopStep1 - startStep1);

  /******************************************************************************************************************************/
  /*                                     Step 2 - Create MST on Clique from step 1                                              */
  /******************************************************************************************************************************/
  startStep2 = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> t1 = g1->PrimMST();

  stopStep2 = std::chrono::high_resolution_clock::now();
  durationStep2 += std::chrono::duration_cast<std::chrono::microseconds>(stopStep2 - startStep2);

  /******************************************************************************************************************************/
  /*                               Step 3 - Replace each edge of MST with shortest path                                         */
  /******************************************************************************************************************************/
  startStep3 = std::chrono::high_resolution_clock::now();

  std::vector<uint32_t> treeNodes;
  std::vector<std::shared_ptr<Edge>> treeEdges;
  uint32_t numberOfTreeNodes = 0;
  uint32_t numberOfTreeEdges = 0;
  // przeiteruj cale adj.list t1 i zamien na tmpShortestPaths
  // Trzeba wyekstrachowac wierzcholki i krawędzie
  // BEZ POWTÓRZEN KRAW I WIERZCH
  for(uint32_t i = 0; i < numberOfSteinerPoints; ++i) {
    for (uint32_t j = 0; j < t1->adjacencyList[i].size(); ++j) {
      uint32_t start = t1->adjacencyList[i].at(j)->start->id;
      uint32_t end = t1->adjacencyList[i].at(j)->end->id;
      int32_t idx = findInEdgeVector(start, end, ShortestPaths);
      if (idx > -1) {
        for (uint32_t k = 0; k < tmpShortestPaths.at(idx)->size(); ++k) {
          // TODO make part below more transparent possibly add addEdgeIfNotAlreadyIn
          // we can do that bc ShortestPaths.at(i) corresponds to tmpShortestPaths[i]
          std::shared_ptr<Edge> e = tmpShortestPaths.at(idx)->at(k);
          int32_t repetitionCheck = findInEdgeVector(e->start->id, e->end->id, treeEdges);
          if (repetitionCheck == -1) {
            treeEdges.push_back(e);
            ++numberOfTreeEdges;

            addNodeIfNotAlreadyIn(e->start->id, treeNodes, numberOfTreeNodes);
            addNodeIfNotAlreadyIn(e->end->id, treeNodes, numberOfTreeNodes);
          }
        }
      } else {
        std::cerr << "Error: couldnt find edge in shortestPath augmention" << std::endl;
        printEdge(t1->adjacencyList[i].at(j));
      }
    }
  }

  std::shared_ptr<Graph> gs(new Graph(treeNodes, treeEdges, numberOfTreeNodes, numberOfTreeEdges, printFlag));

  stopStep3 = std::chrono::high_resolution_clock::now();
  durationStep3 += std::chrono::duration_cast<std::chrono::microseconds>(stopStep3 - startStep3);

  /******************************************************************************************************************************/
  /*                                     Step 4 - Create MST on graph from step 3                                              */
  /******************************************************************************************************************************/
  startStep4 = std::chrono::high_resolution_clock::now();

  std::shared_ptr<Graph> ts = gs->PrimMST();

  stopStep4 = std::chrono::high_resolution_clock::now();
  durationStep4 += std::chrono::duration_cast<std::chrono::microseconds>(stopStep4 - startStep4);

  /******************************************************************************************************************************/
  /*                     Step 5 - Delete all non-terminal leafs untill all leafs remaining are terminals                        */
  /******************************************************************************************************************************/
  startStep5 = std::chrono::high_resolution_clock::now();

  std::vector<std::shared_ptr<Edge>> tsTreeEdges = flattenAdjacencyList(ts->adjacencyList, ts->numberOfNodes);

  uint32_t numberOfSteinerTreeNodes     = numberOfTreeNodes;
  uint32_t tmpNumberOfSteinerTreeNodes  = numberOfTreeNodes;

  uint32_t numberOfSteinerTreeEdges     = tsTreeEdges.size();
  uint32_t tmpNumberOfSteinerTreeEdges  = tsTreeEdges.size();

  startNotNecessary = std::chrono::high_resolution_clock::now();

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfTreeNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(ts, localCopyOfAdjacencyList);

  stopNotNecessary = std::chrono::high_resolution_clock::now();
  durationNotNecessary += std::chrono::duration_cast<std::chrono::microseconds>(stopNotNecessary - startNotNecessary);

  bool foundDanglindBranch = false;
  bool continueSearch = true;
  while (continueSearch) {
    for (uint32_t i = 0; i < numberOfSteinerTreeNodes; ++i) {
      //if non terminal leaf

      //TODO check if it's safe (e.i. if start == localCopyOfAdjacencyList[i])
      if (localCopyOfAdjacencyList[i].size() == 1 && findInUintVector(localCopyOfAdjacencyList[i].at(0)->start->id, terminals) == -1) {
        foundDanglindBranch = true;

        treeNodes.erase(std::remove(treeNodes.begin(), treeNodes.end(), localCopyOfAdjacencyList[i].at(0)->start->id), treeNodes.end());
        --tmpNumberOfSteinerTreeNodes;

        std::shared_ptr<Edge> e = findInEdgeVectorAndReturnValue(
          localCopyOfAdjacencyList[i].at(0)->start,
          localCopyOfAdjacencyList[i].at(0)->end,
          tsTreeEdges);
        if (e != nullptr) {
          tsTreeEdges.erase(std::remove(tsTreeEdges.begin(), tsTreeEdges.end(), e), tsTreeEdges.end());
          --tmpNumberOfSteinerTreeEdges;
          removeEdgeFromAdjacencyList(e, localCopyOfAdjacencyList);
        } else {
          ts->printAdajcencyListFromGraph();
          compareAdajcencyLists(ts->adjacencyList,localCopyOfAdjacencyList, ts->numberOfNodes);
          std::cerr << "localCopyOfAdjacencyList["<<i<<"].at(0)->start " << localCopyOfAdjacencyList[i].at(0)->start->id << std::endl;
          std::cerr << "localCopyOfAdjacencyList["<<i<<"].at(0)->end " << localCopyOfAdjacencyList[i].at(0)->end->id << std::endl;
          printEdgeVector(localCopyOfAdjacencyList[i]);
          std::cerr << "ts->vertices["<<i<<"] " << ts->vertices[i]->id << std::endl;
          std::cerr << "ts->numberOfNodes " << ts->numberOfNodes << std::endl;
          std::cerr << "numberOfTreeNodes " << numberOfTreeNodes << std::endl;
          printEdgeVector(tsTreeEdges);
          printUintVector(*nodeArrayToUint(ts->vertices, ts->numberOfNodes));

          std::cerr << "Edge not found in delete leaves part" << std::endl;
          return dummySharedPointerGraph();
        }
      }
    }
    if (!foundDanglindBranch) {
      continueSearch = false;
      // break;
    } else {
      foundDanglindBranch = false;
      sortAdjacencyList(numberOfSteinerTreeNodes, localCopyOfAdjacencyList);
      numberOfSteinerTreeNodes = tmpNumberOfSteinerTreeNodes;
      numberOfSteinerTreeEdges = tmpNumberOfSteinerTreeEdges;
    }
  }
  resetVisitedStatus();
  std::shared_ptr<Graph> th(new Graph(treeNodes, tsTreeEdges, numberOfSteinerTreeNodes, numberOfSteinerTreeEdges, printFlag));

  stopStep5 = std::chrono::high_resolution_clock::now();
  durationStep5 += std::chrono::duration_cast<std::chrono::microseconds>(stopStep5 - startStep5);

  /***                                             Time for all the prints                                                     ***/
  startNotNecessary = std::chrono::high_resolution_clock::now();

  if (printFlag && numberOfNodes < 1000) {
    std::cout << "KOU" << std::endl;
    for (uint32_t i = 0; i < th->numberOfNodes; ++i) {
      for (uint32_t j = 0; j < th->adjacencyList[i].size(); ++j) {
        std::cout << th->adjacencyList[i].at(j)->start->id << "->" <<   th->adjacencyList[i].at(j)->end->id << "; ";
      }
    }
  }
  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;

  stopNotNecessary = std::chrono::high_resolution_clock::now();
  durationNotNecessary += std::chrono::duration_cast<std::chrono::microseconds>(stopNotNecessary - startNotNecessary);

  timeMeasurements.push_back(durationNotNecessary);
  timeMeasurements.push_back(durationStep1);
  timeMeasurements.push_back(durationStep2);
  timeMeasurements.push_back(durationStep3);
  timeMeasurements.push_back(durationStep4);
  timeMeasurements.push_back(durationStep5);
  return th;
}
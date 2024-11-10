#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <functional>

/*
Lists should have the same number of vectors.
@return returns true if they are different
*/
bool compareAdajcencyLists(std::vector<std::shared_ptr<Edge>>* adjList1, std::vector<std::shared_ptr<Edge>>* adjList2, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    if (adjList1[i].size() != adjList2[i].size()) {
      std::cerr << "Error: Different number of elementsa at "<< i << std::endl;
      for (uint32_t j = 0; j < adjList1[i].size(); ++j)
        std::cout << adjList1[i].at(j)->end->id << ",";
      std::cout << std::endl;
      for (uint32_t j = 0; j < adjList2[i].size(); ++j)
        std::cout << adjList2[i].at(j)->end->id << ",";
      std::cout << std::endl;
      return true;
    } else {
      for (uint32_t j = 0; j < adjList1[i].size(); ++j) {
        if (adjList1[i].at(j)->end->id != adjList2[i].at(j)->end->id && adjList1[i].at(j)->end->id != adjList2[i].at(j)->end->id)
        {
          std::cerr << "Error: Difference in row " << i << " at element "<< j << " of value: "
            << adjList1[i].at(j)->start->id << ","<< adjList1[i].at(j)->end->id << " vs "
            << adjList2[i].at(j)->start->id << ","<< adjList2[i].at(j)->end->id  << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}

void copyAdjacencyListFromGraphWithNewNodeInstances(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList) {
  for (uint32_t i = 0; i < graph->numberOfNodes; ++i)
    for (const std::shared_ptr<Edge>& edge : graph->adjacencyList[i]) {
      int32_t idx = findInArray(edge->start->id, graph->vertices, graph->numberOfNodes);
      if (idx != -1) {
        std::shared_ptr<Node> start(new Node(graph->vertices[idx]->id));
        idx = findInArray(edge->end->id, graph->vertices, graph->numberOfNodes);
        if (idx != -1) {
          std::shared_ptr<Node> end(new Node(graph->vertices[idx]->id));
          copyOfAdjacencyList[i].push_back(std::shared_ptr<Edge>(new Edge(start, edge->weight, end)));
        } else {
          std::cerr << "Error: Couldn't find edge during copyAdjacencyListFromGraphWithNewNodeInstances: " << std::endl; printEdge(edge);
        }
      } else {
        std::cerr << "Error: Couldn't find edge during copyAdjacencyListFromGraphWithNewNodeInstances: " << std::endl; printEdge(edge);
      }
    }
}

/*
returns index from vector which holds value
*/
int32_t findInUintVector(uint32_t value, std::vector<uint32_t> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if (vec.at(i) == value)
      return i;
  return -1;
}

/*
returns index from vector which holds edge
*/
int32_t findInEdgeVector(uint32_t start, uint32_t end, std::vector<std::shared_ptr<Edge>> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if ((vec.at(i)->start->id == start && vec.at(i)->end->id == end) || (vec.at(i)->start->id == end && vec.at(i)->end->id == start))
      return i;
  return -1;
}

/*
returns index from vector which holds edge
*/
std::shared_ptr<Edge> findInEdgeVectorAndReturnValue(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if ((vec.at(i)->start == start && vec.at(i)->end == end) || (vec.at(i)->start == end && vec.at(i)->end == start))
      return vec.at(i);
  return std::shared_ptr<Edge>(nullptr);
}


/*
returns idex from array of Nodes which holds value
*/
int32_t findInArray(uint32_t value, std::shared_ptr<Node>* array, uint32_t arraySize) {
  for (uint32_t i = 0; i < arraySize; ++i)
    if (array[i]->id == value)
      return i;
  return -1;
}

/*
returns idex from array of Nodes which holds value
*/
int32_t findInTmpTree(uint32_t edgeStart, uint32_t edgeEnd, std::vector<PseudoEdge> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if ((vec.at(i).start == edgeStart && vec.at(i).end == edgeEnd) || (vec.at(i).start == edgeEnd && vec.at(i).end == edgeStart))
      return i;
  return -1;
}

/*
finds an edge in adjacencyList and returns pointer to it
*/
std::shared_ptr<Edge> findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<std::shared_ptr<Edge>>* adjacencyList) {

  for(uint32_t i = 0; i < adjacencyList[edgeStart].size(); ++i) {
    if (adjacencyList[edgeStart].at(i) == nullptr) {
      std::cerr << "Error: adjacencyList[" << edgeStart << "].at(" << i << ") returns nullpointer" << std::endl;
      return nullptr;
    }
    if (adjacencyList[edgeStart].at(i)->end->id == edgeEnd) //TODO should be ok, bc we have 2 copies of this edge
      return adjacencyList[edgeStart].at(i);
  }
  std::cerr << "Error: findEdge returns nullpointer for "<< edgeStart << " - " << edgeEnd << std::endl;
  return nullptr;
}

/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void resetVisitedStatusInCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList,  std::shared_ptr<Graph> g) {
  for (uint32_t i = 0; i < g->numberOfNodes; ++i) {
    for (uint32_t j = 0; j < localCopyOfAdjacencyList[i].size(); ++j) {
      if (localCopyOfAdjacencyList[i].at(j)->weight)
        localCopyOfAdjacencyList[i].at(j)->weight = g->adjacencyList[i].at(j)->weight;
      localCopyOfAdjacencyList[i].at(j)->end->visited = false;
      localCopyOfAdjacencyList[i].at(j)->start->visited = false;
    }
  }
}

/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void fullResetCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList,  std::shared_ptr<Graph> g) {
  for (uint32_t i = 0; i < g->numberOfNodes; ++i) {
    for (uint32_t j = 0; j < localCopyOfAdjacencyList[i].size(); ++j) {
      if (localCopyOfAdjacencyList[i].at(j)->weight)
        localCopyOfAdjacencyList[i].at(j)->weight = g->adjacencyList[i].at(j)->weight;
      localCopyOfAdjacencyList[i].at(j)->end->visited = false;
      localCopyOfAdjacencyList[i].at(j)->start->visited = false;
      localCopyOfAdjacencyList[i].at(j)->pred = nullptr;
      localCopyOfAdjacencyList[i].at(j)->succ = nullptr;
    }
  }
}

/*
Finds both instances of edge and updates their weights.
*/
void updateWeight(std::shared_ptr<Edge>& e, uint32_t weight, std::vector<std::shared_ptr<Edge>>*& adjList) {
  e->weight = weight;
  std::shared_ptr<Edge> e2 = findEdge(e->end->id, e->start->id, adjList);
  if (e2 == nullptr) {
    std::cerr << "Error: Couldn't find edge: "; printEdge(e);
  } else {
    e2->weight = weight;
  }
}

/*
Finds both instances of edge and updates their pred.
*/
void updatePred(std::shared_ptr<Edge>& e, std::shared_ptr<Edge> pred, std::vector<std::shared_ptr<Edge>>* adjList) {
  e->pred = pred;
  std::shared_ptr<Edge> e2 = findEdge(e->end->id, e->start->id, adjList);
  if (e2 != nullptr) {
    e2->pred = pred;
  } else {
    std::cerr << "Error: Edge not found, cannot set pred." << std::endl;
  }
}

/*
Finds both instances of edge and creates pred loop.
*/
void updatePredToLoop(std::shared_ptr<Edge>& e, std::vector<std::shared_ptr<Edge>>* adjList) {
  e->pred = e;
  std::shared_ptr<Edge> e2 = findEdge(e->end->id, e->start->id, adjList);
  if (e2 != nullptr) {
    e2->pred = e2;
  } else {
    std::cerr << "Error: Edge not found, cannot set pred." << std::endl;
  }
}

/*
Finds unique numbers in vector and returns number of them
*/
uint32_t tmpPseudoEdgeFindUniqueNumbers(std::vector<PseudoEdge> vec) {
  std::vector<uint32_t> v;
  for (const PseudoEdge& p : vec) {
    v.push_back(p.start);
    v.push_back(p.end);
  }
  std::sort(v.begin(), v.end());
  return std::unique(v.begin(), v.end()) - v.begin();
}

/*
Finds unique numbers in vector and returns vetor of them
*/
std::vector<uint32_t> tmpPseudoEdgeReturnUniqueNumbers(std::vector<PseudoEdge> vec) {
  std::vector<uint32_t> v;
  for (const PseudoEdge& p : vec) {
    v.push_back(p.start);
    v.push_back(p.end);
  }
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

void addTotmptreeEdgesIfNotAlreadyIn(std::vector<PseudoEdge>& tmptreeEdges, std::shared_ptr<Edge> edg) {
  if (findInTmpTree(edg->start->id, edg->end->id, tmptreeEdges) == -1)
    tmptreeEdges.push_back(PseudoEdge(edg->start->id, edg->end->id, edg->weight));
}

/*
Maybe add field to the Edge - isNULL?
*/
std::shared_ptr<Edge> findZeroEdgeInAdjacentTo(uint32_t idx, std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList) {
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[idx].size(); ++i)
    if (localCopyOfAdjacencyList[idx].at(i)->weight == 0)
      return localCopyOfAdjacencyList[idx].at(i);
  std::cerr << "Error: findZeroEdgeInAdjacentTo returns nullptr; searched index: " << idx << std::endl;
  return nullptr;
}

/*
Returns number of edges in clique of |V| = numberOfNodes
*/
uint32_t numberOfEdgesInClique(uint32_t numberOfNodes) {
  uint32_t sum = 0;
  for (uint32_t i = 1; i < numberOfNodes; ++i)
    sum += i;
  return sum;
}

//TODO WROOOOOOOONG WE HAVE TO FIND WHICH ROW HAS start->id AS VALUE!!!!!!!
/*
Removes single edge instance from adjacency list, used as helper function for removeEdgeFromAdjacencyList
*/
void removeSingleEdgeFromAdjacencyList(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>>*& adjList) {
  std::shared_ptr<Edge> edge = findInEdgeVectorAndReturnValue(start, end, adjList[start->id]);
  std::cout << "removing edge: " << std::endl;
  printEdge(edge);
  if (edge != nullptr) {
    adjList[start->id].erase(std::remove(adjList[start->id].begin(), adjList[start->id].end(), edge), adjList[start->id].end());
  } else {
    std::cerr << "Error: Edge not found in removeFromAdjacencyList part" << std::endl;
  }
}

/*
Removes both instances of edge from adjacency list
*/
void removeEdgeFromAdjacencyList(std::shared_ptr<Edge> e, std::vector<std::shared_ptr<Edge>>*& adjList) {
  removeSingleEdgeFromAdjacencyList(e->start, e->end, adjList);
  removeSingleEdgeFromAdjacencyList(e->end, e->start, adjList);
}

/*
Moves non empty vectors to free spaces in the front of the array
*/
void sortAdjacencyList(uint32_t currentAdjListSize, std::vector<std::shared_ptr<Edge>>*& adjList) {
  //TODO ponder if optimizeable
  for (uint32_t i = 0; i < currentAdjListSize; ++i) {
    if (adjList[i].size() == 0) {
      for (uint32_t j = i + 1; j < currentAdjListSize; ++j) {
        if (adjList[j].size() != 0) {
          //TODO ponder if it's safe
          std::swap(adjList[i], adjList[j]);
          break;
        }
      }
    }
  }
}

/*
Adds node to given vector if it is not already in and increases counter
*/
void addNodeIfNotAlreadyIn(uint32_t nodeId, std::vector<uint32_t>& vec, uint32_t& counter) {
  int32_t repetitionCheck = findInUintVector(nodeId, vec);
  if (repetitionCheck == -1) {
    vec.push_back(nodeId);
    ++counter;
  }
}

/*
Adds edge to given vector if it is not already in and increases counter
*/
void addEdgeIfNotAlreadyIn(std::shared_ptr<Edge> edge, std::vector<std::shared_ptr<Edge>>& vec, uint32_t& counter) {
  int32_t repetitionCheck = findInEdgeVector(edge->start->id, edge->end->id, vec);
  if (repetitionCheck == -1) {
    vec.push_back(edge);
    ++counter;
  }
}

// TODO posegreguj te funkcje jakos sensownie
/*
w
*/
void embedEdgeIntoTree(
    std::shared_ptr<Edge> &e,
    std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
    std::vector<std::shared_ptr<Edge>>* &adjacencyList,
    std::vector<PseudoEdge> &tmpTreeEdges) {
  updateWeight(e, 0, localCopyOfAdjacencyList);
  updatePredToLoop(e, localCopyOfAdjacencyList);
  std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
  addTotmptreeEdgesIfNotAlreadyIn(tmpTreeEdges, edg);
}
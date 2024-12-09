#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <functional>
#include <unordered_set>

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
Moves non empty vectors to free spaces in the front of the array
*/
void sortAdjacencyList(uint32_t currentAdjListSize, std::vector<std::shared_ptr<Edge>>*& adjList) {
  for (uint32_t i = 0; i < currentAdjListSize; ++i) {
    if (adjList[i].size() == 0) {
      for (uint32_t j = i + 1; j < currentAdjListSize; ++j) {
        if (adjList[j].size() != 0) {
          std::swap(adjList[i], adjList[j]);
          break;
        }
      }
    }
  }
}

// Hash function for Edge based on node IDs
struct EdgeHash {
  std::size_t operator()(const std::shared_ptr<Edge>& edge) const {
    return std::hash<uint32_t>()(edge->start->id) ^ std::hash<uint32_t>()(edge->end->id);
  }
};

// Equality comparison for Edge, an edge 1 - 2 in this sense in equal to an edge 2 - 1
struct EdgeEqual {
  bool operator()(const std::shared_ptr<Edge>& edge1, const std::shared_ptr<Edge>& edge2) const {
    return (edge1->start->id == edge2->start->id && edge1->end->id == edge2->end->id) ||
           (edge1->start->id == edge2->end->id && edge1->end->id == edge2->start->id);
  }
};

std::vector<std::shared_ptr<Edge>> flattenAdjacencyList(
    std::vector<std::shared_ptr<Edge>>* adjacencyList,
    uint32_t numberOfNodes) {

  std::vector<std::shared_ptr<Edge>> result;
  std::unordered_set<std::shared_ptr<Edge>, EdgeHash, EdgeEqual> seenEdges;

  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::vector<std::shared_ptr<Edge>> adjacencyVector = adjacencyList[i];
    for (const auto& edge : adjacencyVector) {
      // Check if either the edge or its reverse has already been added
      if (seenEdges.find(edge) == seenEdges.end()) {
        result.push_back(edge);
        seenEdges.insert(edge);
      }
    }
  }
  return result;
}

/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void resetVisitedStatusAndWeightInCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList,  std::shared_ptr<Graph> g) {
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
  for (uint32_t i = 0; i < vec.size(); ++i) {
    if ((vec.at(i)->start->id == start->id && vec.at(i)->end->id == end->id) ||
        (vec.at(i)->start->id == end->id && vec.at(i)->end->id == start->id))
      return vec.at(i);
  }
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
  // std::cerr << "Error: findEdge returns nullpointer for "<< edgeStart << " - " << edgeEnd << std::endl;
  return nullptr;
}

/*
  Returns instance of edge of weight equal to 0 adjacenct to node of given index.
*/
std::shared_ptr<Edge> findZeroEdgeInAdjacentTo(uint32_t idx, std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList) {
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[idx].size(); ++i)
    if (localCopyOfAdjacencyList[idx].at(i)->weight == 0)
      return localCopyOfAdjacencyList[idx].at(i);
  std::cerr << "Error: findZeroEdgeInAdjacentTo returns nullptr; searched index: " << idx << std::endl;
  return nullptr;
}

uint32_t findShortestPathAndReturnWeight(
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> shortestPathEdgeVec,
    uint32_t startId,
    uint32_t endId) {

  for (std::shared_ptr<std::vector<std::shared_ptr<Edge>>>& shortestPath : shortestPathEdgeVec) {
    if (shortestPath->size() < 2) {
      if (shortestPath->size() == 1) {
        Edge e = (*(*shortestPath).at(0));
        if ((e.start->id == startId && e.end->id == endId) ||
            (e.end->id == startId && e.start->id == endId))
        {
          return (*(*shortestPath).at(0)).weight;
        } else {
          continue;
        }
      } else {
        std::cerr << "Path length = " << shortestPath->size() << std::endl;
        return std::numeric_limits<uint32_t>::max();
      }
    } else {
      Edge e0 = (*(*shortestPath).at(0));
      Edge e1 = (*(*shortestPath).at(1));

      Edge ePrevToLast = (*(*shortestPath).at(shortestPath->size() - 2));
      Edge eLast       = (*(*shortestPath).at(shortestPath->size() - 1));

      // search for startId if found do nothing, else skip
      if (e0.start->id == startId && !(e1.start->id == startId || e1.end->id == startId)) {
      } else if (e0.end->id == startId && !(e1.start->id == startId || e1.end->id == startId)) {
      } else {
        if (eLast.start->id == startId && !(ePrevToLast.start->id == startId || ePrevToLast.end->id == startId)) {
        } else if (eLast.end->id == startId && !(ePrevToLast.start->id == startId || ePrevToLast.end->id == startId)) {
        } else { // startId not found
          continue;
        }
      }

      //search for endId if found return path weight, else skip
      if (e0.start->id == endId && !(e1.start->id == endId || e1.end->id == endId)) {
        return calculatePathWeight(*shortestPath);
      } else if (e0.end->id == endId && !(e1.start->id == endId || e1.end->id == endId)) {
        return calculatePathWeight(*shortestPath);
      } else {
        if (eLast.start->id == endId && !(ePrevToLast.start->id == endId || ePrevToLast.end->id == endId)) {
          return calculatePathWeight(*shortestPath);
        } else if (eLast.end->id == endId && !(ePrevToLast.start->id == endId || ePrevToLast.end->id == endId)) {
          return calculatePathWeight(*shortestPath);
        } else {
          // endId not found
          // continue bo moze byc wiele z sciezek ze startId
          continue;
        }
      }
    }
  }
  std::cerr << "Path from " << startId << " to " << endId << "not found" << std::endl;
  return std::numeric_limits<uint32_t>::max();
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

/*
  Updates weight of given edge to 0.
  Updates its predecessor to loop.
  Finds this edge in given original adjacency list.
  Adds it to a tmpTree of PseudoEdges.
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

std::shared_ptr<std::vector<uint32_t>> nodeArrayToUint(std::shared_ptr<Node>* vertices, uint32_t numberOfNodes) {
  std::shared_ptr<std::vector<uint32_t>> vertexIdVector = std::make_shared<std::vector<uint32_t>>();
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertexIdVector->push_back(vertices[i]->id);
  return vertexIdVector;
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

uint32_t calculatePathWeight(std::vector<std::shared_ptr<Edge>> path) {
  uint32_t pathWeight = 0;
    for (const std::shared_ptr<Edge>& edge : path)
      pathWeight += edge->weight;
    return pathWeight;
}

// Checks if nth bit is set to 1
bool isNthBitSet(uint32_t num, uint32_t n) {
    return (num & (1 << n)) != 0;
}
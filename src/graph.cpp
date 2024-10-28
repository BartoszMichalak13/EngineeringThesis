#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>
#include <stdint.h>

//TODO check resulting trees with bfs
//TODO QoL graph constructor?
//TODO check ranges (e.g. function returning int32 instead of uint32)
//TODO split into files, maybe for debug functions, helpers etc
//TODO remove treeSize?
//TODO inline everything
//TODO powrzucaj consty
//TODO wszystkie uzycie adj powinny byc bezpieczne (tj bfs na steiner tree nie ma w adj[5] wierzcolka 5 itd)
//TODO nazwy, ładne, czytelne, przenoszą sens
//TOOD TakahashiMatsuyama returns 2 trees (tree is not connected)

//coś z resetem pewnie, użyj emplace moze na priorty que (deprecated?)

Random randomGen;
const uint32_t maxEdgeWeight = 1024;

/*
Used in priorityQueues on Edges
*/
struct EdgeWeightComparatorOnPointers {
  bool operator()(std::shared_ptr<Edge> const& e1, std::shared_ptr<Edge> const& e2) {
    return e1->weight > e2->weight;
  }
};

/*
Used in priorityQueues on Edges
*/
struct EdgeWeightComparator {
  bool operator()(Edge const& e1, Edge const& e2) {
    return e1.weight > e2.weight;
  }
};

/*
Struct simulating behaviour of an Edge. Craeted in rage after first attempts at implementations of TM algorithm
*/
struct PseudoEdge
{
  uint32_t start;
  uint32_t end;
  uint32_t weight;
  PseudoEdge(uint32_t start, uint32_t end) {
    this->start = start;
    this->end = end;
  }
  PseudoEdge(uint32_t start, uint32_t end, uint32_t weight) {
    this->start = start;
    this->end = end;
    this->weight = weight;
  }
};

void printEdgePred(std::shared_ptr<Edge> e) {
  if (e->pred == nullptr)
    std::cout << "e.pred doesn't exist; " << std::endl;
  else
    std::cout << "e.pred->end->id " << e->pred->end->id << "; e.pred->start->id " << e->pred->start->id << "; e.pred->weight " << e->pred->weight << "; " << std::endl;
}

void printEdge(std::shared_ptr<Edge> e) {
  if (e == nullptr) {
      std::cout << "e.end doesn't exist" << std::endl;
    return;
  }
  if (e->end == nullptr) {
    std::cout << "e.end doesn't exist" << std::endl;
    // std::cout << "e.end->id" << e.end->id << std::endl;
    return;
  }
  if (e->start == nullptr) {
    std::cout << "e.start doesn't exist" << std::endl;
    return;
  }

  std::cout << "e.end->id " << e->end->id << "; e.start->id " << e->start->id << "; e.weight " << e->weight << "; ";
  if (e->pred == nullptr)
    std::cout << "e.pred doesn't exist; ";
  else
    std::cout << "e.pred->end->id " << e->pred->end->id << "; e.pred->start->id " << e->pred->start->id << "; e.pred->weight " << e->pred->weight << "; ";
  if (e->succ == nullptr)
    std::cout << "e.succ doesn't exist" << std::endl;
  else
    std::cout << "e.succ->end->id " << e->succ->end->id << "; e.succ->start->id " << e->succ->start->id << "; e.succ->weight" << e->succ->weight << std::endl;
}

//TODO RENAME
/*
returns index from vector which holds value
*/
int32_t findInUintVector(uint32_t value, std::vector<uint32_t> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if (vec.at(i) == value)
      return i;
  return -1;
}

//
/*
returns index from vector which holds edge
*/
int32_t findInEdgeVector(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if ((vec.at(i)->start == start && vec.at(i)->end == end) || (vec.at(i)->start == end && vec.at(i)->end == start))
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

// /*
// returns index from vector which holds edge
// */
// std::shared_ptr<Edge> findInVectorAndReturn(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>> vec) {
//   for (uint32_t i = 0; i < vec.size(); ++i)
//     if ((vec.at(i)->start == start && vec.at(i)->end == end) || (vec.at(i)->start == end && vec.at(i)->end == start))
//       return vec.at(i);
//   return std::shared_ptr<Edge>(nullptr);
// }


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
  for(uint32_t i = 0; i < adjacencyList[edgeStart].size(); ++i)
    if (adjacencyList[edgeStart].at(i)->end->id == edgeEnd)
      return adjacencyList[edgeStart].at(i);
  std::cout << "findEdge returns nullpointer for "<< edgeStart << " - " << edgeEnd << std::endl;
  return nullptr;
}

Graph::Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  this->vertices = new std::shared_ptr<Node>[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    vertices[i] = std::shared_ptr<Node>(new Node(i));
}

/*
Used primarly in approxiamtion algorithms
*/
Graph::Graph(std::vector<uint32_t> nodes, std::vector<std::shared_ptr<Edge>> edges, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  this->vertices = new std::shared_ptr<Node>[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertices[i] = std::shared_ptr<Node>(new Node(nodes.at(i)));
  
  for (uint32_t i = 0; i < numberOfEdges; ++i) {
    int32_t idx1 = findInArray(edges.at(i)->start->id, vertices, numberOfNodes);
    if (idx1 < 0)
      std::cout << "Error node not found in Constructor: " << edges.at(i)->start->id << std::endl;
    int32_t idx2 = findInArray(edges.at(i)->end->id, vertices, numberOfNodes);
    if (idx2 < 0)
      std::cout << "Error node not found in Constructor: " << edges.at(i)->end->id << std::endl;
    std::shared_ptr<Node> start = vertices[idx1];
    std::shared_ptr<Node> end = vertices[idx2];
    std::shared_ptr<Edge> e1(new Edge(start, edges.at(i)->weight, end));
    std::shared_ptr<Edge> e2(new Edge(end, edges.at(i)->weight, start));
    this->adjacencyList[idx1].push_back(e1);
    this->adjacencyList[idx2].push_back(e2);
  }
}

Graph::~Graph() {
  delete[] adjacencyList;
  adjacencyList = nullptr;
  delete[] vertices;
  vertices = nullptr;
}

Graph dummyGraph() {
  return Graph(0,0,0);
}

std::shared_ptr<Graph> dummySharedPointerGraph() {
  return std::shared_ptr<Graph>(new Graph(0,0,0));
}

void Graph::addEdge(uint32_t node1Id, uint32_t node2Id) {
  // We choose them, so that they're smaller than numberOfNodes
  const uint32_t weight = randomGen.generateRandomNumber(1, maxEdgeWeight);
  if(weight == 0)
    std::cout << "ERRRRRRRRRRROR\n";
  this->adjacencyList[node1Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node1Id], weight, vertices[node2Id])));
  this->adjacencyList[node2Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node2Id], weight, vertices[node1Id])));
}

void printAdajcencyList(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node "<< i << std::endl;
    std::cout << "\t";
    for (const std::shared_ptr<Edge>& e: adjList[i]) {
      std::cout << e->end->id << " cost(" << e->weight << ")";
      if (e->weight == 0)
        std::cout << "*** ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
}


void Graph::printAdajcencyListFromGraph() {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node " << vertices[i]->id << std::endl;
    std::cout << "\t";
    for (const std::shared_ptr<Edge>& e: adjacencyList[i]) {
      std::cout << e->end->id << " cost(" << e->weight << ")";
      if (e->weight == 0)
        std::cout << "*** ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
}

/*
Lists should have the same number of vectors.
@return returns true if they are different
*/
bool compareAdajcencyLists(std::vector<std::shared_ptr<Edge>>* adjList1, std::vector<std::shared_ptr<Edge>>* adjList2, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    if (adjList1[i].size() != adjList2[i].size()) {
      std::cout << "Different number of elementsa at "<< i << std::endl;
      for (uint32_t j = 0; j < adjList1[i].size(); ++j)
        std::cout << adjList1[i].at(j)->end->id << ",";
      std::cout << std::endl;
      for (uint32_t j = 0; j < adjList2[i].size(); ++j)
        std::cout << adjList2[i].at(j)->end->id << ",";
      std::cout << std::endl;
      return true;
    } else {
      for (uint32_t j = 0; j < adjList1[i].size(); ++j) {
        if (adjList1[i].at(j)->end != adjList2[i].at(j)->end && adjList1[i].at(j)->end != adjList2[i].at(j)->end)
        {
          std::cout << "Difference in row " << i << " at element "<< j << " of value: "
            << adjList1[i].at(j)->start->id << ","<< adjList1[i].at(j)->end->id << " vs "
            << adjList2[i].at(j)->start->id << ","<< adjList2[i].at(j)->end->id  << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}

/*
Generates vector of idx's of terminals.
*/
std::vector<uint32_t> Graph::generateTerminals(uint32_t numberOfTerminals) {
  std::vector<uint32_t> terminals;
  uint32_t term = randomGen.generateRandomNumber(0, numberOfNodes - 1);
  terminals.push_back(term);
  for (uint32_t i = 1; i < numberOfTerminals; ++i) {
    while (findInUintVector(term, terminals) > -1) {
      term = randomGen.generateRandomNumber(0, numberOfNodes - 1);
    }
    terminals.push_back(term);
  }
  return terminals;
}

void Graph::bfs() {
  std::queue<std::shared_ptr<Node>> toVisit;
  std::shared_ptr<Node> currentNode(vertices[0]);
  toVisit.push(currentNode);
  while (!toVisit.empty()) {
    currentNode = toVisit.front();
    toVisit.pop();
    currentNode->visited = true;
    int32_t idx = findInArray(currentNode->id, vertices, numberOfNodes);
    for (const auto& edge : adjacencyList[idx]) {
      if(!edge->end->visited) {
        toVisit.push(edge->end);
        edge->end->visited = true;
      }
    }
  }
} 

/*
copies adjlist of a graph to given copyOfAdjacencyList
*/
void Graph::copyAdjacencyList(std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    for (const std::shared_ptr<Edge>& edge : adjacencyList[i]) {
      std::shared_ptr<Node> start(vertices[edge->start->id]);
      std::shared_ptr<Node> end(vertices[edge->end->id]);
      copyOfAdjacencyList[i].push_back(std::shared_ptr<Edge>(new Edge(start, edge->weight, end, edge->pred, edge->succ)));
    }
}

/*
copies adjlist of a graph to given copyOfAdjacencyList
*/
void copyAdjacencyListFromGraph(std::shared_ptr<Graph>& graph, std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList) {
  for (uint32_t i = 0; i < graph->numberOfNodes; ++i)
    for (const std::shared_ptr<Edge>& edge : graph->adjacencyList[i]) {
      std::shared_ptr<Node> start(graph->vertices[edge->start->id]);
      std::shared_ptr<Node> end(graph->vertices[edge->end->id]);
      copyOfAdjacencyList[i].push_back(std::shared_ptr<Edge>(new Edge(start, edge->weight, end, edge->pred, edge->succ)));
    }
}

/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void resetCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList, Graph * g) {
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
Finds both instances of edge and updates their weights.
*/
void updateWeight(std::shared_ptr<Edge> e, uint32_t weight, std::vector<std::shared_ptr<Edge>>*& adjList) {
  e->weight = weight;
  std::shared_ptr<Edge> edge = findEdge(e->end->id, e->start->id, adjList);
  if (edge == nullptr) {
    std::cout << "Couldn't find edge: "; printEdge(e);
  } else {
    edge->weight = weight;
  }
}

/*
Finds both instances of edge and updates their pred.
*/
void updatePred(std::shared_ptr<Edge> e, std::shared_ptr<Edge> pred, std::vector<std::shared_ptr<Edge>>* adjList) {
  e->pred = pred;
  std::shared_ptr<Edge> e2 = findEdge(e->end->id, e->start->id, adjList);
  if (e2 != nullptr) {
    e2->pred = pred;
  } else {
    std::cout << "Edge not found, cannot set pred." << std::endl;
  }
}

/*
Finds both instances of edge and creates pred loop.
*/
void updatePredToLoop(std::shared_ptr<Edge> e, std::vector<std::shared_ptr<Edge>>* adjList) {
  e->pred = e;
  std::shared_ptr<Edge> e2 = findEdge(e->end->id, e->start->id, adjList);
  if (e2 != nullptr) {
    e2->pred = e2;
  } else {
    std::cout << "Edge not found, cannot set pred." << std::endl;
  }
}

/*
Searches neighbourhood of node at nodeIndex, adds it to toVisitupdatePred
*/
void searchNeighbours(std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> &toVisit, std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList, uint32_t nodeIndex, std::shared_ptr<Edge> currentEdge) {
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[nodeIndex].size(); ++i) {
    std::shared_ptr<Edge> edge = localCopyOfAdjacencyList[nodeIndex].at(i);
    if (edge->weight) // normal procedure
      updatePred(edge, currentEdge, localCopyOfAdjacencyList);
    else // if it's 0 it's part of a tree
      updatePredToLoop(edge, localCopyOfAdjacencyList);
    updateWeight(edge, (edge->weight + currentEdge->weight), localCopyOfAdjacencyList);
    if (edge == nullptr) {
      std::cout << "Null edge detected in searchNeighbours" << std::endl;
    } else {
      toVisit.push(edge);
    }
  }
}

/*
Prints PseudoEdge
*/
void tmpPseudoEdgePrint(PseudoEdge p) {
  std::cout << p.start << " - " << p.end << std::endl;
}

/*
Prints all PseudoEdges in given vector
*/
void tmpPseudoEdgePrintVec(std::vector<PseudoEdge> vec) {
  for (const PseudoEdge& p : vec) {
    tmpPseudoEdgePrint(p);
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

void addTotmptreeEdgesIfNotAlreadyIn(std::vector<PseudoEdge>& tmptreeEdges, std::shared_ptr<Edge> edg, uint32_t& treeSize) {
  if (findInTmpTree(edg->start->id, edg->end->id, tmptreeEdges) == -1) {
    tmptreeEdges.push_back(PseudoEdge(edg->start->id, edg->end->id, edg->weight));
    ++treeSize;
  }
}

/*
Maybe add field to the Edge - isNULL?
*/
std::shared_ptr<Edge> findZeroEdgeInAdjacentTo(uint32_t uniqueNodes, std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList) {
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[uniqueNodes].size(); ++i)
    if (localCopyOfAdjacencyList[uniqueNodes].at(i)->weight == 0)
      return localCopyOfAdjacencyList[uniqueNodes].at(i);
  std::cout << "ERR findZeroEdgeInAdjacentTo returns nullptr " << std::endl;
  return nullptr;
}

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

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  //init
  resetVisitedStatus();
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyList(localCopyOfAdjacencyList);

  // printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);

  std::vector<PseudoEdge> tmptreeEdges;

  uint32_t treeSize = 0;
  bool foundTerminal = false;
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;

  for (uint32_t i = 0; i < localCopyOfAdjacencyList[terminals.at(0)].size(); ++i)
    toVisit.push(localCopyOfAdjacencyList[terminals.at(0)].at(i));

  vertices[terminals.at(0)]->visited = true;
  terminals.erase(terminals.begin());

  do {
    //TODO is this check necessary?
    if (compareAdajcencyLists(adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
      std::cout << "Comparing Adajcency Lists Failed" << std::endl;
      return dummySharedPointerGraph();
    }

    if (foundTerminal) {
      while(!toVisit.empty())
        toVisit.pop();

      std::vector<uint32_t> uniqueNodes = tmpPseudoEdgeReturnUniqueNumbers(tmptreeEdges);
      for (uint32_t i = 0; i < uniqueNodes.size(); ++i) {
        std::shared_ptr<Edge> e = findZeroEdgeInAdjacentTo(uniqueNodes.at(i), localCopyOfAdjacencyList);
        if (e == nullptr) {
          std::cout << "e doesnt exit aka no neighbour with 0 weight for node "<< uniqueNodes.at(i) << std::endl;
          return dummySharedPointerGraph();
        }
        searchNeighbours(toVisit, localCopyOfAdjacencyList, uniqueNodes.at(i), e);
        e = nullptr;
      }
      foundTerminal = false;
    }

    //main loop
    while(!foundTerminal) { // O(n^2)
      //remoeve all edges that lead to already visited nodes
      if (!toVisit.empty())
        while (toVisit.top() == nullptr || toVisit.top()->end->visited)
          toVisit.pop();

      if (toVisit.empty()) {
        std::cout << "End in loop TakahashiMatsuyama; dummy graph" << std::endl;
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
        updateWeight(e, 0, localCopyOfAdjacencyList);
        updatePredToLoop(e, localCopyOfAdjacencyList);
        std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
        addTotmptreeEdgesIfNotAlreadyIn(tmptreeEdges, edg, treeSize);

        // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
        while (oldPred != nullptr && oldPred->start->id != e->start->id) {
          e = oldPred;
          oldPred = e->pred;
          updateWeight(e, 0, localCopyOfAdjacencyList);
          updatePredToLoop(e, localCopyOfAdjacencyList);
          std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
          addTotmptreeEdgesIfNotAlreadyIn(tmptreeEdges, edg, treeSize);
        } // untill we reach beginning of the path
        resetCopyOfAdjacencyList(localCopyOfAdjacencyList, this);
      } else { // aka terminal not found
        searchNeighbours(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
      }
    }
  } while(!terminals.empty()); // O(k)

  uint64_t totalWeight = 0;
  for (uint32_t i = 0; i < treeSize; ++i)
    totalWeight += tmptreeEdges.at(i).weight;

  if (printFlag && numberOfNodes < 1000) {
    std::cout << "Edges in  tmptreeEdges: " << std::endl;
    for (uint32_t i = 0; i < treeSize; ++i)
      std::cout << tmptreeEdges.at(i).start << "->" <<  tmptreeEdges.at(i).end << "; ";
    std::cout << std::endl;
    std::cout << "Original terminals: ";
    for (uint32_t i = 0; i < originalTerminals.size(); ++i)
      std::cout << originalTerminals.at(i) << ", ";
    std::cout << std::endl;
  }
  std::cout << "TakahashiMatsuyama totalWeight = " << totalWeight << std::endl;
  for (uint32_t i = 0; i < originalTerminals.size(); ++i) {
    bool found = false;
    for (uint32_t j = 0; j < treeSize; ++j) {
      if (originalTerminals.at(i) == tmptreeEdges.at(j).start || originalTerminals.at(i) == tmptreeEdges.at(j).end) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << std::endl;
      std::cout << "missing: " << originalTerminals.at(i) << std::endl;
    }
  }
  std::cout << std::endl;
  std::vector<uint32_t> nodes = tmpPseudoEdgeReturnUniqueNumbers(tmptreeEdges);
  std::vector<std::shared_ptr<Edge>> treeEdges;
  for (uint32_t i = 0; i < tmptreeEdges.size(); ++i) {
    std::shared_ptr<Edge> e = findEdge(tmptreeEdges.at(i).start, tmptreeEdges.at(i).end, adjacencyList);
    if (e == nullptr) {
      std::cout << "missing: " << tmptreeEdges.at(i).start << " - " << tmptreeEdges.at(i).end << std::endl;
    }
    treeEdges.push_back(e);
  }
  std::cout << "tmptreeEdges.size() " << tmptreeEdges.size() << " treeSize " << treeSize << std::endl;
  std::shared_ptr<Graph> steinerTree(new Graph(nodes, treeEdges, nodes.size(), tmptreeEdges.size(), printFlag));

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
  return steinerTree;
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
 
/*
This is Dijkstra algorithm, that stops when we have found node2
*/
std::vector<std::shared_ptr<Edge>> Graph::ShortestPath(uint32_t node1, uint32_t node2) {
  std::vector<std::shared_ptr<Edge>> shortestPath;
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyList(localCopyOfAdjacencyList);
  
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[node1].size(); ++i)
    toVisit.push(localCopyOfAdjacencyList[node1].at(i));
  vertices[node1]->visited = true;

  while (!toVisit.empty()) {
    std::shared_ptr<Edge> e = toVisit.top();
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    if (nextNodeIndex == node2) {
      //TODO can we put that code in while part? (same in Takahashi)

      std::shared_ptr<Edge> oldPred = e->pred;
      std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
      shortestPath.push_back(edg);

      // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
      while (oldPred != nullptr && oldPred->start->id != e->start->id) {
        e = oldPred;
        oldPred = e->pred;
        std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
        shortestPath.push_back(edg);
      } // untill we reach beginning of the path
      return shortestPath;
    } else {
      searchNeighbours(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
    }
  }
  std::cerr << "Could not find path between node1 and node2" << std::endl;
  return std::vector<std::shared_ptr<Edge>>(); //empty vec
}

//TODO WROOOOOOOONG WE HAVE TO FIND WHICH ROW HAS start->id AS VALUE!!!!!!!
/*
Removes single edge instance from adjacency list, used as helper function for removeEdgeFromAdjacencyList
*/
void removeSingeEdgeFromAdjacencyList(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>>*& adjList) {
  std::shared_ptr<Edge> edge = findInEdgeVectorAndReturnValue(start, end, adjList[start->id]);
  if (edge != nullptr){
    adjList[start->id].erase(std::remove(adjList[start->id].begin(), adjList[start->id].end(), edge), adjList[start->id].end());
  } else {
    std::cerr << "Edge not found in removeFromAdjacencyList part" << std::endl;
  }
}

/*
Removes both instances of edge from adjacency list
*/
void removeEdgeFromAdjacencyList(std::shared_ptr<Edge> e, std::vector<std::shared_ptr<Edge>>*& adjList) {
  removeSingeEdgeFromAdjacencyList(e->start, e->end, adjList);
  removeSingeEdgeFromAdjacencyList(e->end, e->start, adjList);
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
  int32_t repetitionCheck = findInEdgeVector(edge->start, edge->end, vec);
  if (repetitionCheck == -1) {
    vec.push_back(edge);
    ++counter;
  }
}

//TODO na spokojnie sprawdz czy all pointery dzialaja tj, czy sa dobrze przepisywane
/*
1. Construct complete graph G_1 where V = terminals (Steiner Points), E = shortest paths between terminals
2. Find mst of G_1
3. Construct graph G_s by replacing each edge from mst with corresponding shortest path
4. Find mst of G_s
5. Remove all branches that lead to leaf which is not Steiner Point (terminal)

In all steps above if there are multiple answers to choose from (e.g. multiple msts) choos aribitrary one

Possible imporvements: what if we create not a clique, but a "normal" graph, we use same method, but we do not contract edges
We may check basic algorithm and my changed version
Most likely (almost obvoiusly) it will only change execution time, not resulting weight of Steiner Tree
*/
std::shared_ptr<Graph> Graph::KouMarkowskyBerman(std::vector<uint32_t> terminals) {

  const uint32_t numberOfSteinerPoints = terminals.size();
  const uint32_t numberOfCliqueEdges = numberOfEdgesInClique(numberOfSteinerPoints);

  std::vector<std::shared_ptr<Edge>> ShortestPaths;
  std::vector<std::shared_ptr<Edge>>* tmpShortestPaths = new std::vector<std::shared_ptr<Edge>>[numberOfCliqueEdges];

  uint32_t currentShortestPathIdx = 0;
  for (uint32_t i = 0; i < numberOfSteinerPoints; ++i) {
    for (uint32_t j = i + 1; j < numberOfSteinerPoints; ++j) {
      tmpShortestPaths[currentShortestPathIdx] = ShortestPath(terminals.at(i), terminals.at(j));
      ++currentShortestPathIdx;
    }
  }

  for (uint32_t i = 0; i < numberOfCliqueEdges; ++i) {
    uint32_t shortestPathWeight = 0;
    for (uint32_t j = 0; j < tmpShortestPaths[i].size(); ++j) {
      shortestPathWeight += tmpShortestPaths[i].at(j)->weight;
    }
    std::shared_ptr<Node> start = tmpShortestPaths[i].at(0)->start;
    std::shared_ptr<Node> end = tmpShortestPaths[i].at(tmpShortestPaths[i].size() - 1)->end;
    ShortestPaths.push_back(std::shared_ptr<Edge>(new Edge(start, shortestPathWeight, end)));
  }

  std::shared_ptr<Graph> g1(new Graph(terminals, ShortestPaths, numberOfSteinerPoints, numberOfCliqueEdges, printFlag));

  std::shared_ptr<Graph> t1 = g1->PrimMST();

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
        for (uint32_t k = 0; k < tmpShortestPaths[idx].size(); ++k) {
          //TODO make part below more transparent possibly add addEdgeIfNotAlreadyIn
          //we can do that bc ShortestPaths.at(i) corresponds to tmpShortestPaths[i]
          std::shared_ptr<Edge> e = tmpShortestPaths[idx].at(k);
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

  std::shared_ptr<Graph> gs(new Graph(treeNodes, treeEdges, numberOfTreeNodes, numberOfTreeEdges, printFlag));

  std::shared_ptr<Graph> ts = gs->PrimMST();

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
  resetVisitedStatus();
  std::shared_ptr<Graph> th(new Graph(treeNodes, treeEdges, numberOfSteinerTreeNodes, numberOfSteinerTreeEdges, printFlag));
  return th;
}

/*
Initialize a tree with a single vertex, chosen arbitrarily from the graph.
Grow the tree by one edge: Of the edges that connect the tree to vertices not yet in the tree, find the minimum-weight edge, and transfer it to the tree.
Repeat step 2 (until all vertices are in the tree).
*/
std::shared_ptr<Graph> Graph::PrimMST() {// do it with priority Queue, maybe my own immplementation?{
  resetVisitedStatus();
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;
  for (std::shared_ptr<Edge>& edge : adjacencyList[0]) //init
    toVisit.push(edge);
  vertices[0]->visited = true;

  std::vector<uint32_t> nodes;
  std::vector<std::shared_ptr<Edge>> tmptreeEdges;

  // std::shared_ptr<Edge> tmptreeEdges[numberOfNodes - 1];
  uint64_t totalWeight = 0;
  for (uint32_t treeSize = 0; treeSize < numberOfNodes - 1; ++treeSize) {
    while (toVisit.top()->end->visited) //remoeve all edges that lead to already visited nodes
      toVisit.pop();

    if (toVisit.size() == 0) {
      std::cout << "End in loop PimMST; dummy graph" << std::endl;
      return dummySharedPointerGraph();
    }
    std::shared_ptr<Edge> e = toVisit.top();
    tmptreeEdges.push_back(e);

    int32_t repetitionCheck = findInUintVector(e->start->id, nodes);
    if (repetitionCheck == -1) {
      nodes.push_back(e->start->id);
    }

    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;

    repetitionCheck = findInUintVector(nextNodeIndex, nodes);
    if (repetitionCheck == -1) {
      nodes.push_back(nextNodeIndex);
    }

    for (std::shared_ptr<Edge>& edge : adjacencyList[nextNodeIndex]) {
      if (!edge->end->visited) {
        toVisit.push(edge);
      }
    }
    totalWeight += e->weight;
  }
  if (printFlag && numberOfNodes < 1000) {
    for (uint32_t i = 0; i < numberOfNodes - 1; ++i)
      std::cout << tmptreeEdges.at(i)->start->id << "->" <<  tmptreeEdges.at(i)->end->id << "; ";
    std::cout << std::endl;
  }
  std::cout << "PimMST totalWeight = " << totalWeight << std::endl;
  resetVisitedStatus();
  return std::shared_ptr<Graph>(new Graph(nodes, tmptreeEdges, numberOfNodes, numberOfNodes - 1, printFlag));
}

void Graph::printData() {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node " << i << " is connected to:\n\t"; 
    for (const auto& edge : adjacencyList[i]) {
      std::cout << static_cast<uint32_t>(edge->end->id) << " cost(" << edge->weight << ") ";
      if(edge->weight == 0)
        std::cout << "ERROR;";
    }
    std::cout << std::endl;
  }
}

void Graph::resetVisitedStatus() {
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertices[i]->visited = false;
}

void Graph::printVisitedStatus() {
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    std::cout << "node " << i << " visited = " << vertices[i]->visited << std::endl;
}

bool Graph::isConnected() {
    std::cout << "HELLO in" << std::endl;

  resetVisitedStatus();
    std::cout << "HELLO in" << std::endl;

  bfs();
      std::cout << "HELLO in" << std::endl;

  bool returnValue = true;
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i]->visited)
      return false;
  std::cout << "HELLO in" << std::endl;

  resetVisitedStatus();
      std::cout << "HELLO in" << std::endl;

  return returnValue;
}

bool Graph::checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id) {
  for (const auto& edge : adjacencyList[node1Id])
    if (edge->end->id == node2Id)
      return true;
  return false;
}

void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool  printFlag)
{
  Graph graph = Graph(numberOfNodes, numberOfEdges, printFlag);
  const uint32_t range = numberOfNodes - 1;

  for (uint32_t i = 0; i < numberOfEdges; ++i)
  {
    bool cannotMakeEdge = true;
    while(cannotMakeEdge) 
    {
      uint32_t node1Id = randomGen.generateRandomNumber(0, range);
      uint32_t node2Id = randomGen.generateRandomNumber(0, range);
      cannotMakeEdge = graph.checkIfEdgeExists(node1Id, node2Id);
      if (!cannotMakeEdge && node1Id != node2Id) //no multi edges, no self-loops
        graph.addEdge(node1Id, node2Id);
    }
  }
  if (printFlag && numberOfNodes < 1000)
  {
    graph.printData();
  }

  if(graph.isConnected()) {
    std::cout << "Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    return;
  }
  std::shared_ptr<Graph> mst = graph.PrimMST(); //dummy for now
  graph.resetVisitedStatus();

  if(mst->isConnected()) { // TODO move code below to isConnected
    std::cout << "mst is connected" << std::endl;
  } else {
    std::cout << "mst is NOT connected" << std::endl;
    return;
  }

  uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
  std::vector<uint32_t> terminals = graph.generateTerminals(numberOfTerminals);

  std::shared_ptr<Graph> steinerTreeTakahashiMatsuyama(graph.TakahashiMatsuyama(terminals));
    std::cout << "HELLO" << std::endl;

  graph.resetVisitedStatus();
  std::cout << "HELLO" << std::endl;

  if(steinerTreeTakahashiMatsuyama->isConnected()) { // TODO move code below to isConnected
    std::cout << "steinerTreeTakahashiMatsuyama is connected" << std::endl;
  } else {
    std::cout << "steinerTreeTakahashiMatsuyama is NOT connected" << std::endl;
    // return;
  }

  std::cout << "HELLO" << std::endl;

  std::shared_ptr<Graph> steinerTreeKouMarkowskyBerman(graph.KouMarkowskyBerman(terminals));
  graph.resetVisitedStatus();

  if(steinerTreeKouMarkowskyBerman->isConnected()) { // TODO move code below to isConnected
    std::cout << "steinerTreeKouMarkowskyBerman is connected" << std::endl;
  } else {
    std::cout << "steinerTreeKouMarkowskyBerman is NOT connected" << std::endl;
    return;
  }

}

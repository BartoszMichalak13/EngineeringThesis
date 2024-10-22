#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>
#include <stdint.h>

//TODO segfault na na koniec; do vec kopiuje sie podwojnie, zle kopiuje wierzcholki z wektora, nie wszystkie, a potem problem ze nie znajduje indexu w tablicy / za duzy
//TODO Graf nie jest polaczony

//coś z resetem pewnie, użyj emplace moze na priorty que

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

/*
returns index from vector which holds value
*/
int32_t findInVector(uint32_t value, std::vector<uint32_t> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if (vec.at(i) == value)
      return i;
  return -1;
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
int32_t findInTmpTree(uint32_t edgeStart, uint32_t edgeEnd, std::vector<PseudoEdge> array) {
  for (uint32_t i = 0; i < array.size(); ++i)
    if ((array.at(i).start == edgeStart && array.at(i).end == edgeEnd) || (array.at(i).start == edgeEnd && array.at(i).end == edgeStart))
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
//Original terminals: x3, x23, x39, x11, x25, x32, x17, x30, x20, x6,
Graph::Graph(std::vector<std::shared_ptr<Edge>> edges, std::vector<uint32_t> nodes, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  std::cout << "numberOfNodes: " << numberOfNodes << std::endl;
  std::cout << "numberOfEdges: " << numberOfEdges << std::endl;
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  this->vertices = new std::shared_ptr<Node>[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i)
  {
    vertices[i] = std::shared_ptr<Node>(new Node(nodes.at(i)));
    std::cout << "vertices: " << nodes.at(i) << " vs " << vertices[i]->id << std::endl;
  }
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
  // adjacencyList = nullptr;
  delete[] vertices;
  // vertices = nullptr;
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
    while (findInVector(term, terminals) > -1) {
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
    for (const auto& edge : adjacencyList[currentNode->id]) {
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
      // Edge copyOfEdge = new Edge(start, edge.weight, end, edge.pred, edge.succ);
      std::shared_ptr<Edge> copyOfEdgePtr(new Edge(start, edge->weight, end, edge->pred, edge->succ));// = &copyOfEdge;
      copyOfAdjacencyList[i].push_back(copyOfEdgePtr);
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
    // std::cout << "added " << std::endl;
    // printEdge(edg);
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
2. in this variation of Dijkstra we search just untill we meet node we search for
3. Terminals should be narrowed down after each found one
KIEDY ZNAJDE WAGI W NOWYM GRAFIE NA 0000000000000000000

usuwam wiekszosc pointerow, gdyz tworze nowe drzewo, a nie chce uszkadzac starego

modify edges, give pointer to pred

@param terminals should be of positive size
*/
std::shared_ptr<Graph> Graph::Dijkstra(std::vector<uint32_t> terminals) {// do it with priority Queue, maybe my own immplementation?

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  //init
  resetVisitedStatus();
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyList(localCopyOfAdjacencyList);

  printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);

  std::vector<PseudoEdge> tmptreeEdges;

  uint32_t treeSize = 0;
  bool foundTerminal = false;
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;

  for (uint32_t i = 0; i < localCopyOfAdjacencyList[terminals.at(0)].size(); ++i)
    toVisit.push(localCopyOfAdjacencyList[terminals.at(0)].at(i));

  vertices[terminals.at(0)]->visited = true;
  terminals.erase(terminals.begin());

  do {
    if (compareAdajcencyLists(adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
      std::cout << "Comparing Adajcency Lists Failed" << std::endl;
      return dummySharedPointerGraph();
    }

    if (foundTerminal) {
      while(!toVisit.empty())
        toVisit.pop();

      uint32_t currentNumberOfNodes = tmpPseudoEdgeFindUniqueNumbers(tmptreeEdges);
      std::vector<uint32_t> uniqueNodes = tmpPseudoEdgeReturnUniqueNumbers(tmptreeEdges);
      std::cout << "currentNumberOfNodes " << currentNumberOfNodes << std::endl;
      std::cout << "debug tsize " << treeSize << std::endl;
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
        std::cout << "End in loop Dijkstra; dummy graph" << std::endl;
        return dummySharedPointerGraph();
      }
      std::shared_ptr<Edge> e = toVisit.top();
      toVisit.pop();
      e->end->visited = true;
      uint32_t nextNodeIndex = e->end->id;

      // check if we have found terminal
      int32_t idx = findInVector(nextNodeIndex, terminals);
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
  std::cout << "Dijkstra totalWeight = " << totalWeight << std::endl;

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
  std::shared_ptr<Graph> steinerTree(new Graph(treeEdges, nodes, nodes.size(), tmptreeEdges.size(), printFlag));

  std::cout << "dupa 321" << std::endl;
  for (uint32_t i = 0; i < steinerTree->numberOfNodes; ++i)
    std::cout << "v["<< i <<"]= "  << steinerTree->vertices[i]->id << std::endl;
  std::cout << "dupa 321" << std::endl;

  //TODO NIE DELETOWAC ALE TRZEBA
  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;

  for (uint32_t i = 0; i < steinerTree->numberOfNodes; ++i)
    std::cout << "v["<< i <<"]= "  << steinerTree->vertices[i]->id << std::endl;
  std::cout << "dupa 321" << std::endl;
  return steinerTree;
}

/*
Initialize a tree with a single vertex, chosen arbitrarily from the graph.
Grow the tree by one edge: Of the edges that connect the tree to vertices not yet in the tree, find the minimum-weight edge, and transfer it to the tree.
Repeat step 2 (until all vertices are in the tree).
*/
Graph Graph::PrimMST() // do it with priority Queue, maybe my own immplementation?
{
  resetVisitedStatus();
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;
  for (std::shared_ptr<Edge>& edge : adjacencyList[0]) //init
    toVisit.push(edge);
  vertices[0]->visited = true;
  
  std::shared_ptr<Edge> tmptreeEdges[numberOfNodes - 1];
  uint64_t totalWeight = 0;
  for (uint32_t treeSize = 0; treeSize < numberOfNodes - 1; ++treeSize)
  {
    while (toVisit.top()->end->visited) //remoeve all edges that lead to already visited nodes
      toVisit.pop();

    if (toVisit.size() == 0) 
    {
      std::cout << "End in loop PimMST; dummy graph" << std::endl;
      return dummyGraph();
    }
    std::shared_ptr<Edge> e = toVisit.top();
    tmptreeEdges[treeSize] = e;
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    for (std::shared_ptr<Edge>& edge : adjacencyList[nextNodeIndex])
    {
      if (!edge->end->visited)
      {
        toVisit.push(edge);
      }
    }
    totalWeight += e->weight;
  }
  if (printFlag && numberOfNodes < 1000)
  {
    for (uint32_t i = 0; i < numberOfNodes - 1; ++i)
      std::cout << tmptreeEdges[i]->start->id << "->" <<  tmptreeEdges[i]->end->id << "; ";
    std::cout << std::endl;
  }
  std::cout << "PimMST totalWeight = " << totalWeight << std::endl;

  std::cout << "Dummy graph primMST" << std::endl;
  return dummyGraph();
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

bool Graph::isConnected() {
  bfs();
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i]->visited)
      return false;
  return true;
}

void Graph::resetVisitedStatus() {
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertices[i]->visited = false;
}

void Graph::printVisitedStatus() {
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    std::cout << "node " << i << " visited = " << vertices[i]->visited << std::endl;
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
    std::cout << "dupa" << std::endl;

  if(graph.isConnected()) {
    std::cout << "Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    //TMP SOLUTION
    return;
  }
      std::cout << "dupa" << std::endl;

  graph.resetVisitedStatus();
  Graph mst = graph.PrimMST();
  graph.resetVisitedStatus();
        std::cout << "d2314upa" << std::endl;

  mst.printAdajcencyListFromGraph();

  //TODO: maybe make them not random?
  //TMP SOLUTION
  uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
  std::cout << "dupa" << std::endl;
  std::vector<uint32_t> terminals = graph.generateTerminals(numberOfTerminals);
  std::cout << "dupa" << std::endl;

  // Graph steinerTree = graph.Dijkstra(terminals);
  std::shared_ptr<Graph> steinerTree(graph.Dijkstra(terminals));
  std::cout << "dupa" << std::endl;
  std::cout << "steinerTree.numberOfNodes " << steinerTree->numberOfNodes << std::endl;
  std::cout << "steinerTree.numberOfEdges " << steinerTree->numberOfEdges << std::endl;
  for (uint32_t i = 0; i < steinerTree->numberOfNodes; ++i)
    std::cout << "v["<< i <<"]= "  << steinerTree->vertices[i]->id << std::endl;
  std::cout << "dupa" << std::endl;
  std::cout << "steinerTree.adjacencyList[0].size() " << steinerTree->adjacencyList[0].size() << std::endl;
  std::cout << "graph.adjacencyList[0].size() " << graph.adjacencyList[0].size() << std::endl;
  graph.printAdajcencyListFromGraph();
    std::cout << "dupa" << std::endl;
  steinerTree->printAdajcencyListFromGraph();
  std::cout << "dupa" << std::endl;

  for (uint32_t i = 0; i < steinerTree->numberOfNodes; ++i)
    for (uint32_t j = 0; j < steinerTree->adjacencyList[i].size(); ++j) {
      printEdge(steinerTree->adjacencyList[i].at(j));
      printEdgePred(steinerTree->adjacencyList[i].at(j));
    }
  std::cout << "dupa" << std::endl;
}



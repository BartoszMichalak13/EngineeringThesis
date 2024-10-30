#include "graph.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "random.hpp"

#include <iostream>
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
//TOOD TakahashiMatsuyama didnt find all nodes
//TOOD TakahashiMatsuyama psuje liczbe wierzcholkow w orginalnym grafie
//TODO potencjalnie mozna sprawdzac czy acykliczny

//coś z resetem pewnie, użyj emplace moze na priorty que (deprecated?)

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

Graph Graph::dummyGraph() {
  return Graph(0,0,0);
}

std::shared_ptr<Graph> Graph::dummySharedPointerGraph() {
  return std::shared_ptr<Graph>(new Graph(0,0,0));
}

void Graph::addEdge(uint32_t node1Id, uint32_t node2Id) {
  // We choose them, so that they're smaller than numberOfNodes

  const uint32_t weight = Random::getInstance()->generateRandomNumber(1, maxEdgeWeight);
  if(weight == 0)
    std::cout << "ERRRRRRRRRRROR\n";
  this->adjacencyList[node1Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node1Id], weight, vertices[node2Id])));
  this->adjacencyList[node2Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node2Id], weight, vertices[node1Id])));
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
Generates vector of idx's of terminals.
*/
std::vector<uint32_t> Graph::generateTerminals(uint32_t numberOfTerminals) {
  std::vector<uint32_t> terminals;
  uint32_t term = Random::getInstance()->generateRandomNumber(0, numberOfNodes - 1);
  terminals.push_back(term);
  for (uint32_t i = 1; i < numberOfTerminals; ++i) {
    while (findInUintVector(term, terminals) > -1) {
      term = Random::getInstance()->generateRandomNumber(0, numberOfNodes - 1);
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

// /*
// copies adjlist of a graph to given copyOfAdjacencyList
// */
// void Graph::copyAdjacencyList(std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList) {
//   for (uint32_t i = 0; i < numberOfNodes; ++i)
//     for (const std::shared_ptr<Edge>& edge : adjacencyList[i]) {
//       std::shared_ptr<Node> start(vertices[edge->start->id]);
//       std::shared_ptr<Node> end(vertices[edge->end->id]);
//       copyOfAdjacencyList[i].push_back(std::shared_ptr<Edge>(new Edge(start, edge->weight, end, edge->pred, edge->succ)));
//     }
// }


//TODO ideas: moze glowne adjacency list ma nie wyczyszczone pred?
/*
Searches neighbourhood of node at nodeIndex, adds it to toVisitupdatePred
*/
void Graph::searchNeighbours(
      std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> &toVisit,
      std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
      uint32_t nodeIndex,
      std::shared_ptr<Edge> currentEdge) {

  int32_t idx = findInArray(nodeIndex, vertices, numberOfNodes);
  if (idx == -1) {
    std::cerr << "Error: Node " << nodeIndex << " not found in searchNeighbours"  << std::endl;
    return;
  }

  for (uint32_t i = 0; i < localCopyOfAdjacencyList[idx].size(); ++i) {
    std::shared_ptr<Edge> edge = localCopyOfAdjacencyList[idx].at(i);
    // if they are different edges
    if (!((edge->start->id == currentEdge->start->id  && edge->end->id == currentEdge->end->id) ||
        (edge->start->id == currentEdge->end->id    && edge->end->id == currentEdge->start->id))) {
      if (edge->weight) // normal procedure
        updatePred(edge, currentEdge, localCopyOfAdjacencyList); //TODO check if works ok
      else // if it's 0 it's part of a tree
        updatePredToLoop(edge, localCopyOfAdjacencyList); //TODO: I think it should be already done

      if (edge->weight) // if positive update, else does not change, it's tree  part KEY PART
        updateWeight(edge, (edge->weight + currentEdge->weight), localCopyOfAdjacencyList);
      if (edge == nullptr) {
        std::cerr << "Error: Null edge detected in searchNeighbours" << std::endl;
      } else {
        toVisit.push(edge);
      }
    } else {
      // std::cerr << "Error: SearchNeighbours currentEdge is the same as edge" << std::endl;
      // printEdge(edge);
      // std::cerr << "currentEdge" << std::endl;
      // printEdge(currentEdge);
    }
  }
}
//TODO in kou in shortespath part add edges from adjList, not local adj list
 
/*
This is Dijkstra algorithm, that stops when we have found node2
*/
std::shared_ptr<std::vector<std::shared_ptr<Edge>>> Graph::ShortestPath(uint32_t node1, uint32_t node2) {
  std::shared_ptr<Graph> self = shared_from_this();

  // std::shared_ptr<std::vector<std::shared_ptr<Edge>>> shortestPath(std::vector<std::shared_ptr<Edge>>());
  std::shared_ptr<std::vector<std::shared_ptr<Edge>>> shortestPath = std::make_shared<std::vector<std::shared_ptr<Edge>>>();

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraph(self, localCopyOfAdjacencyList);

  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;
  int32_t idx = findInArray(node1, vertices, numberOfNodes);
  if (idx == -1) {
    std::cerr << "Error: Node " << node1 << " not found in ShortestPath"  << std::endl;
    return std::shared_ptr<std::vector<std::shared_ptr<Edge>>>();
  }
  for (uint32_t i = 0; i < localCopyOfAdjacencyList[idx].size(); ++i)
    toVisit.push(localCopyOfAdjacencyList[idx].at(i));
  vertices[idx]->visited = true;

  while (!toVisit.empty()) {
    std::shared_ptr<Edge> e = toVisit.top();
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    if (nextNodeIndex == node2) {
      //TODO can we put that code in while part? (same in Takahashi)
      std::shared_ptr<Edge> oldPred = e->pred;
      std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
      (*shortestPath).push_back(edg);

      // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
      // while (oldPred != nullptr && oldPred->start->id != e->start->id) {
      while (oldPred != nullptr && e->start->id != node1) {
        e = oldPred;
        oldPred = e->pred;
        std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
        printEdge(e);
        printEdgePred(e);
        (*shortestPath).push_back(edg);
      } // untill we reach beginning of the path
      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      return shortestPath;
    } else {
      searchNeighbours(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
    }
    // std::cout << "HELLO del bef" << std::endl;
    // while ((toVisit.top() == nullptr || toVisit.top()->end->visited) && !toVisit.empty()) {//TODO is it always the end?
    //   if (toVisit.top() == nullptr) {
    //    std::cout << "HELLO del in nullptr" << std::endl;

    //   } else {
    //      std::cout << "HELLO del in ";
    //     printEdge(toVisit.top());
    //   }
    //   toVisit.pop();
    // }
    // std::cout << "HELLO del aft" << std::endl;
  }
  std::cerr << "Error: Could not find path between node1 and node2" << std::endl;
  delete[] localCopyOfAdjacencyList;
  localCopyOfAdjacencyList = nullptr;
  return std::shared_ptr<std::vector<std::shared_ptr<Edge>>>(); //empty vec
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

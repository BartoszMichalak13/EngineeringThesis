#include "graph.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include "random.hpp"

#include <iostream>
#include <functional>
#include <stdint.h>

//TODO QoL graph constructor?
//TODO check ranges (e.g. function returning int32 instead of uint32)
//TODO split into files, maybe for debug functions, helpers etc
//TODO inline everything
//TODO powrzucaj consty
//TODO wszystkie uzycie adj powinny byc bezpieczne (tj bfs na steiner tree nie ma w adj[5] wierzcolka 5 itd)
//TODO nazwy, ładne, czytelne, przenoszą sens
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

std::shared_ptr<Graph> dummySharedPointerGraph() {
  return std::shared_ptr<Graph>(new Graph(0,0,0));
}

/*
Used only for the first graph, else we may want to change vertices[nodeId] to idx = find(NodeId); vertices[idx]
*/
void Graph::addEdgeWithRandomWeight(const uint32_t node1Id, const uint32_t node2Id) {
  // We choose them, so that they're smaller than numberOfNodes
  const uint32_t weight = Random::getInstance()->generateRandomNumber(1, maxEdgeWeight);
  this->adjacencyList[node1Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node1Id], weight, vertices[node2Id])));
  this->adjacencyList[node2Id].push_back(std::shared_ptr<Edge>(new Edge(vertices[node2Id], weight, vertices[node1Id])));
}

/*
Used only for the first graph, else we may want to change vertices[nodeId] to idx = find(NodeId); vertices[idx]
*/
void Graph::addEdge(const uint32_t node1Id, const uint32_t node2Id, const uint32_t weight) {
  // We choose them, so that they're smaller than numberOfNodes
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

// void Graph::bfs() {
//   std::queue<std::shared_ptr<Node>> toVisit;
//   std::shared_ptr<Node> currentNode(vertices[0]);
//   toVisit.push(currentNode);
//   while (!toVisit.empty()) {
//     currentNode = toVisit.front();
//     toVisit.pop();
//     currentNode->visited = true;
//     int32_t idx = findInArray(currentNode->id, vertices, numberOfNodes);
//     for (const auto& edge : adjacencyList[idx]) {
//       if(!edge->end->visited) {
//         toVisit.push(edge->end);
//         edge->end->visited = true;
//       }
//     }
//   }
// }

bool Graph::bfs() {
  bool isAcyclic = true;
  std::queue<std::shared_ptr<Node>> toVisit;
  std::shared_ptr<Node> currentNode(vertices[0]);
  toVisit.push(currentNode);
  while (!toVisit.empty()) {
    currentNode = toVisit.front();
    toVisit.pop();
    if (currentNode->visited) {
      isAcyclic = false;
    }
    currentNode->visited = true;
    int32_t idx = findInArray(currentNode->id, vertices, numberOfNodes);
    for (const auto& edge : adjacencyList[idx]) {
      if (!edge->end->visited) {
        toVisit.push(edge->end);
      }
    }
  }
  return isAcyclic;
} 

/*
Calculates total cost of all graph edges
*/
uint32_t Graph::graphTotalCost() {
  uint32_t cost = 0;
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    for (const auto& edge : adjacencyList[i]) {
      cost += edge->weight;
      if(edge->weight == 0)
        std::cerr << "ERROR;" << std::endl;
    }
  }
  return cost/2; //we counted both edges 0-1 and 1-0 
}

//TODO clear this mess
//TODO zle value dla sciezki miedzy x i x (petla)
/*
Searches neighbourhood of node at nodeIndex, adds it to toVisitupdatePred
*/
void Graph::searchNeighboursV2(
      std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> &toVisit,
      std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
      uint32_t nodeIndex,
      std::shared_ptr<Edge> currentEdge)
{

  int32_t idx = findInArray(nodeIndex, vertices, numberOfNodes);
  if (idx == -1) {
    std::cerr << "Error: Node " << nodeIndex << " not found in searchNeighboursV2"  << std::endl;
    return;
  }

  for (uint32_t i = 0; i < localCopyOfAdjacencyList[idx].size(); ++i)
  {
    std::shared_ptr<Edge> edge = localCopyOfAdjacencyList[idx].at(i);
    // if they are different edges
    if (!((edge->start->id == currentEdge->start->id  && edge->end->id == currentEdge->end->id) ||
        (edge->start->id == currentEdge->end->id    && edge->end->id == currentEdge->start->id)))
    {
      if(!edge->end->visited) {
        if (edge->weight)
        {// normal procedure
          if (edge->pred == nullptr)
          { // if pred != nullptr it means the edge was seen before, and may be a part of another path
            updatePred(edge, currentEdge, localCopyOfAdjacencyList);
          }
          std::shared_ptr<Edge> e = findEdge(edge->start->id, edge->end->id, adjacencyList);
          if (e == nullptr)
          {
            std::cerr << "Error: edge e is nullptr in adjacencyList in searchNeighboursV2" << std::endl;
            return;
          }
          //it is crucial to take original edge weight + current, as if the edge'd have been seen before the weight'd been much greater
          updateWeight(edge, (e->weight + currentEdge->weight), localCopyOfAdjacencyList);

        }
        else
        {// if it's 0 it's part of a tree
          updatePredToLoop(edge, localCopyOfAdjacencyList); //TODO: I think it should be already done
        }


        // if (edge->weight)
        // {// if positive update, else does not change, it's tree  part KEY PART
        //   std::shared_ptr<Edge> e = findEdge(edge->start->id, edge->end->id, adjacencyList);
        //   if (e == nullptr)
        //   {
        //     std::cerr << "Error: edge e is nullptr in adjacencyList in searchNeighboursV2" << std::endl;
        //     return;
        //   }
        //   //it is crucial to take original edge weight + current, as if the edge'd have been seen before the weight'd been much greater
        //   updateWeight(edge, (e->weight + currentEdge->weight), localCopyOfAdjacencyList);
        // }


        // if (edge == nullptr)
        // {
        //   std::cerr << "Error: Null edge detected in searchNeighboursv2" << std::endl;
        // }
        // else
        // {
          toVisit.push(edge);
        // }
      }
    }
  }
}

/*
All pairs shortest paths
vector of vectors of Edges

TODO does it always work as intended?
*/
// std::shared_ptr<std::vector<std::shared_ptr<Edge>>>* Graph::AllPairsShortestPath(std::vector<uint32_t> terminals) {
std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> Graph::AllPairsShortestPath(std::vector<uint32_t> terminals) {

  // get graph instance
  std::shared_ptr<Graph> self = shared_from_this();

  // make copy of terminals
  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));


  // create copy of adjacency list (to make changes in this copy and not in original instance of adjacency list)
  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(self, localCopyOfAdjacencyList);

  // create priority queue
  std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> toVisit;

  // create structure for all pairs of shortest paths
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> allPairsOfShortestPaths(0);

  for ( uint32_t currentStartingNodeIndex = 0;
        currentStartingNodeIndex < originalTerminals.size();
        ++currentStartingNodeIndex) {

    // create shortest path structure
    std::shared_ptr<std::vector<std::shared_ptr<Edge>>> shortestPath = std::make_shared<std::vector<std::shared_ptr<Edge>>>();

    /*
      This is the main node from which Dijsktra algorithm starts
      It's always 0 as we truncate terminals on the go
    */
    uint32_t startingNode = terminals.at(0);
    terminals.erase(terminals.begin() + 0);

    // find instance of startingNode in vertices
    int32_t startingNodeIdx = findInArray(startingNode, vertices, numberOfNodes);
    if (startingNodeIdx == -1) {
      std::cerr << "Error: Node " << startingNode << " not found in ShortestPath"  << std::endl;
      return std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>>();
    }

    //add all neighbours of startingNode
    for (uint32_t i = 0; i < localCopyOfAdjacencyList[startingNodeIdx].size(); ++i)
      toVisit.push(localCopyOfAdjacencyList[startingNodeIdx].at(i));
    vertices[startingNodeIdx]->visited = true; //TODO does it always work?

    while (!toVisit.empty() && !terminals.empty()) {
      std::shared_ptr<Edge> e = toVisit.top();
      toVisit.pop();
      e->end->visited = true;
      uint32_t nextNodeIndex = e->end->id;

      int32_t idx = findInUintVector(nextNodeIndex, terminals);
      if (idx > -1) {
        std::shared_ptr<Edge> copyOfe = e;
        terminals.erase(terminals.begin() + idx);
        // std::cout << "adding path from " << startingNode << " to " << nextNodeIndex << std::endl;

        //TODO can we put that code in while part? (same in Takahashi)
        std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
        shortestPath->push_back(edg);

        std::shared_ptr<Edge> oldPred = e->pred;
        // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
        // while (oldPred != nullptr && e->start->id != startingNode) {
        while (oldPred != nullptr && (e->start->id != startingNode && e->end->id != startingNode)) {
          e = oldPred;
          oldPred = e->pred;
          std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
          shortestPath->push_back(edg);
        } // untill we reach beginning of the path
        allPairsOfShortestPaths.push_back(shortestPath);

        // reset shortest path
        shortestPath = std::make_shared<std::vector<std::shared_ptr<Edge>>>();

        // reset e to prev state to use it in searchNeighboursV2
        e = copyOfe;
      }
      searchNeighboursV2(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
    }
    if (!terminals.empty()) {
      std::cerr << "Error: Could not find path between " << startingNode << " and ";
      printNodeVector(terminals);

      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      return std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>>();
    } else {

      // reset terminals to their begining state - nodes already checked
      for (uint32_t i = currentStartingNodeIndex + 1; i < originalTerminals.size(); ++i)
        terminals.push_back(originalTerminals.at(i));

      // make it clear for next Dijkstra iteration
      resetVisitedStatus();

      // reset priority queue
      while(!toVisit.empty())
        toVisit.pop();

      // reset localCopyOfAdjacencyList NOTE if this doesn't work, put new localCopyOfAdjacencyList in for
      fullResetCopyOfAdjacencyList(localCopyOfAdjacencyList, self);
      if (printFlag) {
        std::cout << "terminals.size " << terminals.size() << std::endl;
      }
    }

  }
  return allPairsOfShortestPaths;
}

//TODO centralny punkt majacy polaczenie do all innych w grafie, sprawdz jak to sie zachowuje
 
/*
This is Dijkstra algorithm, that stops when we have found node2

Note
*/
std::shared_ptr<std::vector<std::shared_ptr<Edge>>> Graph::ShortestPath(uint32_t node1, uint32_t node2) {
  std::shared_ptr<Graph> self = shared_from_this();
  std::shared_ptr<std::vector<std::shared_ptr<Edge>>> shortestPath = std::make_shared<std::vector<std::shared_ptr<Edge>>>();

  std::vector<std::shared_ptr<Edge>>* localCopyOfAdjacencyList = new std::vector<std::shared_ptr<Edge>>[numberOfNodes];
  copyAdjacencyListFromGraphWithNewNodeInstances(self, localCopyOfAdjacencyList);

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
      while (oldPred != nullptr && e->start->id != node1) {
        e = oldPred;
        oldPred = e->pred;
        std::shared_ptr<Edge> edg = findEdge(e->start->id, e->end->id, adjacencyList);
        (*shortestPath).push_back(edg);
      } // untill we reach beginning of the path
      delete[] localCopyOfAdjacencyList;
      localCopyOfAdjacencyList = nullptr;
      resetVisitedStatus();
      return shortestPath;
    } else {
      searchNeighboursV2(toVisit, localCopyOfAdjacencyList, nextNodeIndex, e);
    }
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

    int32_t idx = findInArray(nextNodeIndex, vertices, numberOfNodes);
    if (idx == -1) {
      std::cerr << "Error: no such a node in verticies: " << nextNodeIndex << std::endl;
      return dummySharedPointerGraph();
    }
    for (std::shared_ptr<Edge>& edge : adjacencyList[idx]) {
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
  if (printFlag) {
    std::cout << "PimMST totalWeight = " << totalWeight << std::endl;
  }
  resetVisitedStatus();
  return std::shared_ptr<Graph>(new Graph(nodes, tmptreeEdges, numberOfNodes, numberOfNodes - 1, printFlag));
}

void Graph::printData() {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node " << vertices[i]->id << " is connected to:\n\t"; 
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

// bool Graph::isConnected() {
//   resetVisitedStatus();
//   bool isAcyclic = bfs();
//   bool returnValue = true;
//   for (uint32_t i = 0; i < numberOfNodes; ++i)
//     if (!vertices[i]->visited)
//       return false;
//   resetVisitedStatus();
//   return returnValue;
// }

/*
Checks if graph is acyclic and connected
.first is true when graph is acyclic
.second is true when graph is connected
*/
std::pair<bool,bool> Graph::isTree() {
  resetVisitedStatus();
  bool isAcyclic = bfs();
  bool isConnected = true;
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i]->visited)
      isConnected = false;
  resetVisitedStatus();
  return std::pair<bool,bool>(isAcyclic, isConnected);
}

bool Graph::checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id) {
  for (const auto& edge : adjacencyList[node1Id])
    if (edge->end->id == node2Id)
      return true;
  return false;
}

//TODO works only on starting graph as adj[src] makes assumption that edge.start == src
std::vector<std::vector<uint32_t>> Graph::toAdjacencyMatrix() {

  // deafault fill with inf
  std::vector<std::vector<uint32_t>> adjacencyMatrix(numberOfNodes, std::vector<uint32_t>(numberOfNodes, std::numeric_limits<uint32_t>::max()));

  // Fill the diagonal with 0s (distance from each vertex to itself)
  for (uint32_t i = 0; i < numberOfNodes; i++) {
    adjacencyMatrix[i][i] = 0;
  }

  // Copy edges from adjacency list to adjacency matrix
  for (uint32_t src = 0; src < numberOfNodes; src++) {
    for (const auto& edge : adjacencyList[src]) {
      // uint32_t start = edge->start->id;
      uint32_t end = edge->end->id;
      uint32_t weight = edge->weight;
      adjacencyMatrix[src][end] = weight;
      adjacencyMatrix[end][src] = weight; // Undirected graph
    }
  }

  return adjacencyMatrix;
}


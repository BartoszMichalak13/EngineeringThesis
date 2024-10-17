#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>
#include <stdint.h>

Random randomGen;
const uint32_t maxEdgeWeight = 1024;

Graph::Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag)
{
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->vertices = new Node[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
  {
    vertices[i] = Node(i);
  }
  this->adjacencyList = new std::vector<Edge>[numberOfNodes];
}

Graph::~Graph()
{
  delete[] adjacencyList;
  delete[] vertices;
}

void Graph::addEdge(uint32_t node1Id, uint32_t node2Id)
{
  // We choose them, so that they're smaller than numberOfNodes
  const uint32_t weight = randomGen.generateRandomNumber(maxEdgeWeight);
  this->adjacencyList[node1Id].push_back(Edge(&vertices[node2Id], weight, &vertices[node1Id]));
  this->adjacencyList[node2Id].push_back(Edge(&vertices[node1Id], weight, &vertices[node2Id]));
}

void Graph::bfs() 
{
  std::queue<Node*> toVisit;
  Node * currentNode = &vertices[0];
  toVisit.push(currentNode);
  while (!toVisit.empty())
  {
    currentNode = toVisit.front();
    toVisit.pop();
    currentNode->visited = true;
    for (const auto& edge : adjacencyList[currentNode->id])
    {
      if(!edge.end->visited)
      {
        toVisit.push(edge.end);
        edge.end->visited = true;
      }
      //printVisitedStatus();
    }
  }
} 

//Q: is really dijkstra? Q2: is it necessery?

// this is an structure which implements the
// operator overloading - it used to be not a pointer
struct EdgeWeightComparator {
    bool operator()(Edge* const& e1, Edge* const& e2)
    {
        return e1->weight > e2->weight;
    }
};
// void Graph::Dijkstra() 
// {
//   // std::queue<Node*> toVisit;
//   std::priority_queue<Edge, std::vector<Edge>, EdgeWeightComparator> toVisit;
//   Edge currentEdge = Edge(&vertices[0], 0, &vertices[0]); //Pseudo edge used just as "starter"
//   toVisit.push(currentEdge);
//   // // Track total cost
//   // uint32_t totalCost = 0;
//   while (!toVisit.empty())
//   {
//     // std::cout << "We went from node " << currentEdge.end->id; // this untrue in the sense that it may show wrong starting dest
//     currentEdge = toVisit.top();
//     toVisit.pop();
//     currentEdge.end->visited = true;
//     // totalCost += currentEdge.weight;
//     // // Print the connection and the weight
//     // std::cout << " to node " << currentEdge.end->id 
//     //           << " with cost " << currentEdge.weight << "\n";
//     for (const auto& edge : adjacencyList[currentEdge.end->id])
//     {
//       if(!edge.end->visited)
//       {
//         toVisit.push(edge);
//         edge.end->visited = true;
//       }
//       //printVisitedStatus();
//     }
//   }
//   // std::cout << "Total cost of the path: " << totalCost << "\n";
// } 

/*
Initialize a tree with a single vertex, chosen arbitrarily from the graph.
Grow the tree by one edge: Of the edges that connect the tree to vertices not yet in the tree, find the minimum-weight edge, and transfer it to the tree.
Repeat step 2 (until all vertices are in the tree).
*/
Graph Graph::PrimMST() // do it with priority Queue, maybe my own immplementation?
{
  resetVisitedStatus();
  std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparator> toVisit;
  uint32_t treeSize = 0;
  for (Edge& edge : adjacencyList[0]) //init 
    toVisit.push(&edge);
  vertices[0].visited = true;
  
  Edge * treeEdges[numberOfNodes - 1];
  uint64_t totalWeight = 0;
  for (uint32_t treeSize = 0; treeSize < numberOfNodes - 1; ++treeSize)
  {
    while (toVisit.top()->end->visited) //remoeve all edges that lead to already visited nodes
      toVisit.pop();

    if (toVisit.size() == 0) 
    {
      std::cout << "End in loop PimMST; dummy graph" << std::endl;
      return Graph(0,0,printFlag);
    }
    Edge * e = toVisit.top();
    treeEdges[treeSize] = e;
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    for (Edge& edge : adjacencyList[nextNodeIndex])
    {
      if (!edge.end->visited)
      {
        toVisit.push(&edge);
      }
    }
    totalWeight += e->weight;
  }
  if (printFlag && numberOfNodes < 1000)
  {
    for (uint32_t i = 0; i < numberOfNodes - 1; ++i)
      std::cout << treeEdges[i]->start->id << "->" <<  treeEdges[i]->end->id << "; ";
    std::cout << std::endl;
  }
  std::cout << "PimMST totalWeight = " << totalWeight << std::endl;

  std::cout << "Dummy graph primMST" << std::endl;
  return Graph(0,0,printFlag);
}

// Graph Graph::PrimMST() // do it with priority Queue, maybe my own immplementation?
// {
//   resetVisitedStatus();
//   std::priority_queue<Edge, std::vector<Edge>, EdgeWeightComparator> toVisit;

//   uint32_t treeSize = 0;
//   Node * treeNodes[numberOfNodes];
//   Edge * treeEdges[numberOfNodes - 1];
//   treeNodes[treeSize] = &vertices[treeSize];

//   uint32_t nextNodeIndex = treeSize;
//   Edge * nextEdge;
//   uint64_t totalWeight = 0;
//   uint32_t minValue = UINT32_MAX;
//   Edge * previousMinEdge;// = &Edge(&Node(),0,&Node()); //dummy 

//   while (treeSize < numberOfNodes - 1)
//   {
//     for (Edge& edge : adjacencyList[nextNodeIndex])
//     {
//       if (edge.weight < minValue && !edge.end->visited)
//       {
//         // previousMinEdge = minValue;
//         minValue = edge.weight;

//         nextNodeIndex = edge.end->id;
//         nextEdge = &edge;
//       }
//     }
//     vertices[nextNodeIndex].visited = true;
//     treeEdges[treeSize] = nextEdge;
//     totalWeight += minValue;
//     ++treeSize;
//   }

//   for (uint32_t i = 0; i < numberOfNodes; ++i)
//     std::cout << treeEdges[i]->start->id << "->" <<  treeEdges[i]->end->id << "; ";
//   std::cout << std::endl;
//   std::cout << "totalWeight = " << totalWeight << std::endl;
// }

void Graph::printData()
{
  for (int i = 0; i < numberOfNodes; ++i)
  {
    std::cout << "node " << i << " is connected to:\n\t"; 
    for (const auto& edge : adjacencyList[i])
    {
      std::cout << static_cast<uint32_t>(edge.end->id) << " cost(" << edge.weight << ") "; 
    }
    std::cout << std::endl;
  }
}

bool Graph::isConnected()
{
  bfs();
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i].visited)
      return false;
  return true;
}

void Graph::resetVisitedStatus()
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertices[i].visited = false;
}

void Graph::printVisitedStatus()
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    std::cout << "node " << i << " visited = " << vertices[i].visited << std::endl;
}

bool Graph::checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id) 
{
  for (const auto& edge : adjacencyList[node1Id])
    if (edge.end->id == node2Id)
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
      uint32_t node1Id = randomGen.generateRandomNumber(range);
      uint32_t node2Id = randomGen.generateRandomNumber(range);
      cannotMakeEdge = graph.checkIfEdgeExists(node1Id, node2Id);
      if (!cannotMakeEdge && node1Id != node2Id) //no multi edges, no self-loops
        graph.addEdge(node1Id, node2Id);
    }
  }
  if (printFlag && numberOfNodes < 1000)
  {
    graph.printData();
  }
  if(graph.isConnected())
    std::cout << "Graph is connected" << std::endl;
  else
    std::cout << "Graph is NOT connected" << std::endl;
  graph.resetVisitedStatus();
  Graph mst = graph.PrimMST();
  graph.resetVisitedStatus();


}



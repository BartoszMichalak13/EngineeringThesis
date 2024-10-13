#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>

Random randomGen;
const uint32_t maxEdgeWeight = 1024;

Graph::Graph(uint32_t numberOfNodes, uint32_t numberOfEdges)
{
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
  this->adjacencyList[node1Id].push_back(Edge(&vertices[node2Id], weight));
  this->adjacencyList[node2Id].push_back(Edge(&vertices[node1Id], weight));
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

// this is an structure which implements the
// operator overloading
struct EdgeWeightComparator {
    bool operator()(Edge const& e1, Edge const& e2)
    {
        // return "true" if "p1" is ordered
        // before "p2", for example:
        return e1.weight < e1.weight;
    }
};
void Graph::Dijkstra() 
{
  // std::queue<Node*> toVisit;
  std::priority_queue<Edge, std::vector<Edge>, EdgeWeightComparator> toVisit;
  Edge currentEdge = Edge(&vertices[0], 0);
  toVisit.push(currentEdge);

  // Track total cost
  uint32_t totalCost = 0;

  while (!toVisit.empty())
  {
    std::cout << "We went from node " << currentEdge.end->id; // this untrue in the sense that it may show wrong starting dest
    currentEdge = toVisit.top();
    toVisit.pop();
    currentEdge.end->visited = true;
    totalCost += currentEdge.weight;
    // Print the connection and the weight
    std::cout << " to node " << currentEdge.end->id 
              << " with cost " << currentEdge.weight << "\n";
    for (const auto& edge : adjacencyList[currentEdge.end->id])
    {
      if(!edge.end->visited)
      {
        toVisit.push(edge);
        edge.end->visited = true;
      }
      //printVisitedStatus();
    }
  }
  std::cout << "Total cost of the path: " << totalCost << "\n";
} 

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

void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density)
{
  Graph graph = Graph(numberOfNodes, numberOfEdges);
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

  graph.printData();
  if(graph.isConnected())
    std::cout << "Graph is connected" << std::endl;
  else
    std::cout << "Graph is NOT connected" << std::endl;
  graph.resetVisitedStatus();
  graph.Dijkstra();
}



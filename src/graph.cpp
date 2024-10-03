#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
Graph::Graph(uint16_t numberOfNodes, uint16_t numberOfEdges)
{
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->vertices = new Node[numberOfNodes];
  for (uint16_t i = 0; i < numberOfNodes; ++i) 
  {
    vertices[i] = Node(i);
  }
  this->adjacencyList = new std::vector<Node*>[numberOfNodes];
}

Graph::~Graph()
{
  delete[] adjacencyList;
  delete[] vertices;
}

void Graph::addEdge(uint16_t node1Id, uint16_t node2Id)
{
  // We choose them, so that they're smaller than numberOfNodes
  this->adjacencyList[node1Id].push_back(&vertices[node2Id]);
  this->adjacencyList[node2Id].push_back(&vertices[node1Id]);
}

void Graph::printData()
{
  for (int i = 0; i < numberOfNodes; ++i)
  {
    std::cout << "node " << i << " is connected to:\n\t"; 
    for (const auto& node : adjacencyList[i])
    {
      std::cout << static_cast<uint16_t>(node->id) << " "; 
    }
    std::cout << std::endl;
  }
}

bool Graph::isConnected()
{
  for (uint16_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i].visited)
      return false;
  return true;
}

void Graph::resetVisitedStatus()
{
  for (uint16_t i = 0; i < numberOfNodes; ++i)
    vertices[i].visited = false;
}

void Graph::printVisitedStatus()
{
  for (uint16_t i = 0; i < numberOfNodes; ++i)
    std::cout << "node " << i << " visited = " << vertices[i].visited << std::endl;
}


// Returns if graph is connected
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
    for (const auto& node : adjacencyList[currentNode->id])
    {
      if(!node->visited)
      {
        toVisit.push(node);
        node->visited = true;
      }
      //printVisitedStatus();
    }
  }
} 

bool Graph::checkIfEdgeExists(uint16_t node1Id, uint16_t node2Id) 
{
  for (const auto& node : adjacencyList[node1Id])
    if (node->id == node2Id)
      return true;
  return false;
}

void generateGraph(uint16_t numberOfNodes, uint16_t numberOfEdges, float density)
{
  Graph graph = Graph(numberOfNodes, numberOfEdges);
  Random randomGen;
  uint16_t range = numberOfNodes - 1;

  for (uint16_t i = 0; i < numberOfEdges; ++i)
  {
    bool cannotMakeEdge = true;
    while(cannotMakeEdge) 
    {
      uint16_t node1Id = randomGen.generateRandomNumber(range);
      uint16_t node2Id = randomGen.generateRandomNumber(range);
      cannotMakeEdge = graph.checkIfEdgeExists(node1Id, node2Id);
      if (!cannotMakeEdge && node1Id != node2Id) //no multi edges, no self-loops
        graph.addEdge(node1Id, node2Id);
    }
  }
  graph.bfs();
  if(graph.isConnected())
    std::cout << "Graph is connected" << std::endl;
  else
    std::cout << "Graph is NOT connected" << std::endl;
  graph.printData();

}



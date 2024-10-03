#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <memory>
#include "node.hpp"

class Graph {
  private:
    uint16_t numberOfNodes;
    Node * vertices;
    uint16_t numberOfEdges;
    std::vector<Node*>* adjacencyList;

  public:
    Graph(uint16_t numberOfNodes, uint16_t numberOfEdges);
    ~Graph();
    void addEdge(uint16_t node1Id, uint16_t node2Id);
    void printData();
    bool checkIfEdgeExists(uint16_t node1Id, uint16_t node2Id);
    void bfs();
    bool isConnected();
    void printVisitedStatus();
    void resetVisitedStatus();
};

void generateGraph(uint16_t numberOfNodes, uint16_t numberOfEdges, float density);

#endif

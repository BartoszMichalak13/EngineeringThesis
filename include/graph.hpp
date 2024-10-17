#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <memory>
#include "node.hpp"

class Graph {
  private:
    uint32_t numberOfNodes;
    Node * vertices;
    uint32_t numberOfEdges;
    std::vector<Edge>* adjacencyList;
    void bfs();
    bool printFlag;
  public:
    Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    ~Graph();
    void addEdge(uint32_t node1Id, uint32_t node2Id);
    void printData();
    bool checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id);
    bool isConnected();
    // void Dijkstra(); deprecated
    Graph PrimMST(); //if it is suppose to return graph we either create childClass named tree or add another constructor for list of edges and nodes to populate the new graph
    // Graph KruskalMST(); deprecated
    void printVisitedStatus();
    void resetVisitedStatus();
};

void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool printFlag);

#endif

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <memory>
#include "node.hpp"

class Graph {
  private:

  public:
    uint32_t numberOfNodes;
    std::shared_ptr<Node>* vertices;
    uint32_t numberOfEdges;
    std::vector<std::shared_ptr<Edge>>* adjacencyList;
    bool printFlag;
    void bfs();
    void printAdajcencyListFromGraph();
    //prev private above

    Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    Graph(std::vector<std::shared_ptr<Edge>> edges, std::vector<uint32_t> nodes, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    ~Graph();

    void addEdge(uint32_t node1Id, uint32_t node2Id);
    void printData();
    bool checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id);
    void copyAdjacencyList(std::vector<std::shared_ptr<Edge>>*&  copyOfAdjacencyList);
    std::vector<uint32_t> generateTerminals(uint32_t numberOfTerminals);

    bool isConnected();
    Graph Dijkstra(std::vector<uint32_t> terminals);
    Graph PrimMST(); //if it is suppose to return graph we either create childClass named tree or add another constructor for list of edges and nodes to populate the new graph
    void printVisitedStatus();
    void resetVisitedStatus();
};

void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool printFlag);

#endif

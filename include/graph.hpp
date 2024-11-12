#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <memory>
#include <queue>
#include <utility>
#include <set>
#include "node.hpp"

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


class Graph : public std::enable_shared_from_this<Graph> { // TODO doczytaj o  std::enable_shared_from_this
  private:

  public:
    // Random randomGen;

    uint32_t numberOfNodes;
    std::shared_ptr<Node>* vertices;
    uint32_t numberOfEdges;
    std::vector<std::shared_ptr<Edge>>* adjacencyList;
    bool printFlag;



    bool bfs();
    void printAdajcencyListFromGraph();
    //prev private above


    Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    Graph(std::vector<uint32_t> nodes, std::vector<std::shared_ptr<Edge>> edges, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    std::shared_ptr<Graph> dummySharedPointerGraph();
    Graph dummyGraph();

    ~Graph();



    void addEdge(uint32_t node1Id, uint32_t node2Id);
    void printData();
    bool checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id);
    uint32_t graphTotalCost();
    std::vector<uint32_t> generateTerminals(uint32_t numberOfTerminals);
    std::shared_ptr<std::vector<std::shared_ptr<Edge>>> ShortestPath(uint32_t node1, uint32_t node2);
    void searchNeighboursV2(
        std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> &toVisit,
        std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
        uint32_t nodeIndex,
        std::shared_ptr<Edge> currentEdge);

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> AllPairsShortestPath(std::vector<uint32_t> terminals);
    std::pair<bool,bool> isTree();
    std::shared_ptr<Graph> TakahashiMatsuyama(std::vector<uint32_t> terminals);
    std::shared_ptr<Graph> KouMarkowskyBerman(std::vector<uint32_t> terminals);
    std::shared_ptr<Graph> PrimMST(); //if it is suppose to return graph we either create childClass named tree or add another constructor for list of edges and nodes to populate the new graph


    uint32_t DreyfusWagner(std::vector<uint32_t> terminals);
    uint32_t calculateSteiner(
        std::vector<uint32_t> C,
        // std::vector<std::vector<uint32_t>> graph,
        std::vector<std::set<uint32_t>> allSubsets,
        uint32_t q);


    std::vector<std::vector<uint32_t>> toAdjacencyMatrix();
    void printMatrix(const std::vector<std::vector<uint32_t>>& matrix);


    void printVisitedStatus();
    void resetVisitedStatus();
};

// void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool printFlag);

#endif

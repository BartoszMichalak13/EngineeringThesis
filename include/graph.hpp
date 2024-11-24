#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <memory>
#include <queue>
#include <utility>
#include <set>
#include <chrono>

#include "GraphStructures.hpp"

const uint32_t maxEdgeWeight = 1024;

/*
Used in priorityQueues on Edges
*/
struct EdgeWeightComparatorOnPointers {
  bool operator()(std::shared_ptr<Edge> const& e1, std::shared_ptr<Edge> const& e2) {
    return e1->weight > e2->weight;
  }
};

// TODO change it to edge
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
  public:
    uint32_t numberOfNodes;
    uint32_t numberOfEdges;

    std::shared_ptr<Node>* vertices;
    std::vector<std::shared_ptr<Edge>>* adjacencyList;

    // Debugs and prints
    bool printFlag;
    void printAdajcencyListFromGraph();
    void printData();

    // Functions that check various things
    bool bfs();
    bool checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id);
    std::pair<bool,bool> isTree();

    //Constructors destructors, dummy instances
    Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    Graph(std::vector<uint32_t> nodes, std::vector<std::shared_ptr<Edge>> edges, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag);
    std::shared_ptr<Graph> dummySharedPointerGraph();
    Graph dummyGraph();
    ~Graph();

    // Generators, adders, etc.
    void addEdgeWithRandomWeight(const uint32_t node1Id, const uint32_t node2Id);
    void addEdge(const uint32_t node1Id, const uint32_t node2Id, const uint32_t weight);
    std::vector<uint32_t> generateTerminals(uint32_t numberOfTerminals);
    void printVisitedStatus();
    void printMatrix(const std::vector<std::vector<uint32_t>>& matrix);

    // Main algorithms
    std::shared_ptr<Graph> TakahashiMatsuyama(
      std::vector<uint32_t> terminals,
      std::vector<std::chrono::microseconds> &timeMeasurements);
    std::shared_ptr<Graph> KouMarkowskyBerman(
      std::vector<uint32_t> terminals,
      std::vector<std::chrono::microseconds> &timeMeasurements);
    std::shared_ptr<Graph> PrimMST();
    uint32_t DreyfusWagner(std::vector<uint32_t> terminals);

    // Support algorihtms and functions
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> AllPairsShortestPath(std::vector<uint32_t> terminals);
    uint32_t calculateSteiner(
        std::vector<uint32_t> C,
        std::vector<std::vector<uint32_t>> adjMatrix,
        std::vector<std::set<uint32_t>> allSubsets,
        uint32_t q);
    uint32_t graphTotalCost();
    std::shared_ptr<std::vector<std::shared_ptr<Edge>>> ShortestPath(uint32_t node1, uint32_t node2);
    void searchNeighboursV2(
        std::priority_queue<std::shared_ptr<Edge>, std::vector<std::shared_ptr<Edge>>, EdgeWeightComparatorOnPointers> &toVisit,
        std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
        uint32_t nodeIndex,
        std::shared_ptr<Edge> currentEdge);
    std::vector<std::vector<uint32_t>> toAdjacencyMatrix();
    void resetVisitedStatus();
    /*
    Removes single edge instance from adjacency list, used as helper function for removeEdgeFromAdjacencyList
    */
    void removeSingleEdgeFromAdjacencyList(
        std::shared_ptr<Node> start,
        std::shared_ptr<Node> end,
        std::vector<std::shared_ptr<Edge>>*& adjList);
    /*
    Removes both instances of edge from adjacency list
    */
    void removeEdgeFromAdjacencyList(std::shared_ptr<Edge> e, std::vector<std::shared_ptr<Edge>>*& adjList);

};

// Various functions that are connected with graph but don't belong to it
std::shared_ptr<Graph> dummySharedPointerGraph();
std::vector<uint32_t> generateTerminals(uint32_t numberOfNodes, uint32_t numberOfTerminals);

#endif

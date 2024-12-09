#ifndef DEBUG_PRINTS_HPP
#define DEBUG_PRINTS_HPP
#include <iostream>
#include <vector>
#include <memory>
#include "graph.hpp"

/*
Prints PseudoEdge.
*/
void tmpPseudoEdgePrint(PseudoEdge p);
/*
Prints all PseudoEdges in given vector.
*/
void tmpPseudoEdgePrintVec(std::vector<PseudoEdge> vec);

/*
Prints start, weight, end of an edge and all the same for it's predecessor (TODO and succ depr)
*/
void printEdge(std::shared_ptr<Edge> e);
/*
Prints start, weight, end of an edge predecessor
*/
void printEdgePred(std::shared_ptr<Edge> e);

/*
Call printEdge on all edges in vector
*/
void printEdgeVector(std::vector<std::shared_ptr<Edge>> vec);
/*
Prints all uints in vector
*/
void printUintVector(std::vector<uint32_t> vec);


/*
Prints adajcency list with costs. Prints *** when 0 is seen for better visibility
*/
void printAdajcencyList(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes);
/*
Prints edge predicessor for each edge in adjacency list
*/
void printAdajcencyListPred(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes);

/*
Used when dealing with vec of vec matrix instance.
Prints values or MAX for std::numeric_limits<uint32_t>::max().
INFO printed values have variable length (are NOT filled with 0s to given size)
*/
void printMatrix(const std::vector<std::vector<uint32_t>>& matrix);

#endif

#ifndef DEBUG_PRINTS_HPP
#define DEBUG_PRINTS_HPP
#include <iostream>

void tmpPseudoEdgePrint(PseudoEdge p);
void tmpPseudoEdgePrintVec(std::vector<PseudoEdge> vec);
void printAdajcencyList(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes);
void printEdgePred(std::shared_ptr<Edge> e);
void printEdge(std::shared_ptr<Edge> e);
void printEdgeVector(std::vector<std::shared_ptr<Edge>> vec);
void printNodeVector(std::vector<uint32_t> vec);


#endif

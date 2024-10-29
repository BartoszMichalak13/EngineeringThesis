#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <iostream>

/*
Prints PseudoEdge
*/
void tmpPseudoEdgePrint(PseudoEdge p) {
  std::cout << p.start << " - " << p.end << std::endl;
}

/*
Prints all PseudoEdges in given vector
*/
void tmpPseudoEdgePrintVec(std::vector<PseudoEdge> vec) {
  for (const PseudoEdge& p : vec) {
    tmpPseudoEdgePrint(p);
  }
}

void printAdajcencyList(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node "<< i << std::endl;
    std::cout << "\t";
    for (const std::shared_ptr<Edge>& e: adjList[i]) {
      std::cout << e->end->id << " cost(" << e->weight << ")";
      if (e->weight == 0)
        std::cout << "*** ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
}

void printEdgePred(std::shared_ptr<Edge> e) {
  if (e->pred == nullptr)
    std::cout << "e.pred doesn't exist; " << std::endl;
  else
    std::cout << "e.pred->end->id " << e->pred->end->id << "; e.pred->start->id " << e->pred->start->id << "; e.pred->weight " << e->pred->weight << "; " << std::endl;
}

void printEdge(std::shared_ptr<Edge> e) {
  if (e == nullptr) {
      std::cout << "e.end doesn't exist" << std::endl;
    return;
  }
  if (e->end == nullptr) {
    std::cout << "e.end doesn't exist" << std::endl;
    return;
  }
  if (e->start == nullptr) {
    std::cout << "e.start doesn't exist" << std::endl;
    return;
  }

  std::cout << "e.start->id " << e->start->id << "; e.end->id " << e->end->id << "; e.weight " << e->weight << "; ";
  if (e->pred == nullptr)
    std::cout << "e.pred doesn't exist; ";
  else
    std::cout <<  "e.pred->start->id " << e->pred->start->id << "; e.pred->end->id " << e->pred->end->id << "; e.pred->weight " << e->pred->weight << "; ";
  if (e->succ == nullptr)
    std::cout << "e.succ doesn't exist" << std::endl;
  else
    std::cout  << "e.succ->start->id " << e->succ->start->id << "; e.succ->end->id " << e->succ->end->id << "; e.succ->weight" << e->succ->weight << std::endl;
}


void printEdgeVector(std::vector<std::shared_ptr<Edge>> vec) {
  std::cout << "print vector" << std::endl;
  for (uint32_t i = 0; i < vec.size(); ++i)
    printEdge(vec.at(i));
}

void printNodeVector(std::vector<uint32_t> vec) {
  std::cout << "print nodes" << std::endl;
  for (uint32_t i = 0; i < vec.size(); ++i)
    std::cout << vec.at(i) << ", ";
  std::cout << std::endl;
}

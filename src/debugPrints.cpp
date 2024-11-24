#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"
#include <iostream>


void tmpPseudoEdgePrint(PseudoEdge p) {
  std::cout << p.start << " - " << p.end << "(" << p.weight << ")" << std::endl;
}

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

void printEdge(std::shared_ptr<Edge> e) {
  if (e == nullptr) {
      std::cout << "e doesn't exist" << std::endl;
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

void printEdgePred(std::shared_ptr<Edge> e) {
  if (e->pred == nullptr)
    std::cout << "e.pred doesn't exist; " << std::endl;
  else
    std::cout << "e.pred->end->id " << e->pred->end->id << "; e.pred->start->id " << e->pred->start->id << "; e.pred->weight " << e->pred->weight << "; " << std::endl;
}

void printEdgeVector(std::vector<std::shared_ptr<Edge>> vec) {
  std::cout << "print vector" << std::endl;
  for (uint32_t i = 0; i < vec.size(); ++i)
    printEdge(vec.at(i));
}

void printUintVector(std::vector<uint32_t> vec) {
  std::cout << "print nodes" << std::endl;
  for (uint32_t i = 0; i < vec.size(); ++i)
    std::cout << "el("<<i<<"):" << vec.at(i) << ", ";
  std::cout << std::endl;
}

void printAdajcencyListPred(std::vector<std::shared_ptr<Edge>>* adjList, uint32_t numberOfNodes) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node "<< i << std::endl;
    std::cout << "\t";
    for (const std::shared_ptr<Edge>& e: adjList[i]) {
      printEdgePred(e);
      std::cout << "edge " << std::endl;
    }
    std::cout << std::endl;
  }
}

void printMatrix(const std::vector<std::vector<uint32_t>>& matrix) {
  for (const auto& row : matrix) {
    for (uint32_t val : row) {
      if (val == std::numeric_limits<uint32_t>::max())
        std::cout << "MAX ";
      else
        std::cout << val << " ";
    }
    std::cout << std::endl;
  }
}
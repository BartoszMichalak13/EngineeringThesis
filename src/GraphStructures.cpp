#include "GraphStructures.hpp"

std::shared_ptr<Edge> constructSelfLoopInitEdge(std::shared_ptr<Node> v) {
  return std::shared_ptr<Edge>(new Edge(v,0,v));
};

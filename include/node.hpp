#ifndef NODE_HPP
#define NODE_HPP

#include <cstdint>

struct Node
{
  uint16_t id;
  bool visited;
  
  Node() : id(0), visited(false) {}
  Node(uint16_t id) : id(id), visited(false) {}
  ~Node() {}

  operator int() const {
    return static_cast<int>(id);
  }
  operator uint16_t() const {
    return id;
  }
};

#endif

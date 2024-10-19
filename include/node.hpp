#ifndef NODE_HPP
#define NODE_HPP

#include <cstdint>

//here we use class instead of struct as we need to derive special nodeclass for adjacency list
class Node {
  public:
    uint32_t id;
    bool visited;
    
    Node() :            id(0),  visited(false) {}
    Node(uint32_t id) : id(id), visited(false) {}
    ~Node() {}

    operator int() const {
      return static_cast<int>(id);
    }
    operator uint32_t() const {
      return id;
    }
};

class Edge{
  public:
    Node * start;
    Node * end;
    uint32_t weight;
    Edge * pred;
    Edge * succ;
    //id?

    Edge(bool dummy)                                                          : start(0),     end(0),   weight(!dummy),  pred(this),    succ(nullptr) {}
    Edge()                                                                    : start(0),     end(0),   weight(0),       pred(nullptr), succ(nullptr) {}
    Edge(Node * end, uint32_t weight, Node * start)                           : start(start), end(end), weight(weight),  pred(nullptr), succ(nullptr) {}
    Edge(Node * end, uint32_t weight, Node * start, Edge * pred)              : start(start), end(end), weight(weight),  pred(pred),    succ(nullptr) {}
    Edge(Node * end, uint32_t weight, Node * start, Edge * pred, Edge * succ) : start(start), end(end), weight(weight),  pred(pred),    succ(succ)    {}

    //operations?
};

#endif

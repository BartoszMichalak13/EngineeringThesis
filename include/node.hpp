#ifndef NODE_HPP
#define NODE_HPP

#include <cstdint>

//here we use class instead of struct as we need to derive special nodeclass for adjacency list
class Node 
{
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

// class WeightedNode : public Node
// {
//   public:
//     uint32_t weight;
//     WeightedNode() :                              Node(),   weight(0) {}
//     WeightedNode(uint32_t id) :                   Node(id), weight(0) {}
//     WeightedNode(uint32_t id, uint32_t weight) :  Node(id), weight(weight) {}
//     ~WeightedNode() {}
// };

class Edge
{
  public:
    Node * start; // it is distinc from end (uint instead of node) because it is not as commonly used
    Node * end;
    uint32_t weight;
    //id?

    // Edge(Node * end) :                                    end(end), weight(0),      start(end) {}
    // Edge(Node * end, uint32_t weight) :                   end(end), weight(weight), start(end) {}
    Edge(Node * end, uint32_t weight,  Node * start) :   end(end), weight(weight), start(start) {}

    //operations?
};

#endif

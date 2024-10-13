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
    Node * end;
    uint32_t weight;
    //id?

    Edge(Node * node) :                   end(node), weight(0) {}
    Edge(Node * node, uint32_t weight) :  end(node), weight(weight) {}

    //operations?
};

#endif

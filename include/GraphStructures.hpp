#ifndef NODE_HPP
#define NODE_HPP

#include <cstdint>
#include <memory>

//here we use class instead of struct as we need to derive special nodeclass for adjacency list
class Node {
  public:
    uint32_t id;
    bool visited;
    
    Node()            : id(0),  visited(false) {}
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
    std::shared_ptr<Node> start;
    uint32_t weight;
    std::shared_ptr<Node> end;

    // Predecessor of this edge
    std::shared_ptr<Edge> pred;

    //TODO unused delete
    std::shared_ptr<Edge> succ;

    Edge(bool dummy)                                                                                                                      : start(0),     weight(!dummy), end(0),   pred(this),    succ(nullptr) {}
    Edge()                                                                                                                                : start(0),     weight(0),      end(0),   pred(nullptr), succ(nullptr) {}
    Edge(std::shared_ptr<Node> start, uint32_t weight, std::shared_ptr<Node> end)                                                         : start(start), weight(weight), end(end), pred(nullptr), succ(nullptr) {}
    Edge(std::shared_ptr<Node> start, uint32_t weight, std::shared_ptr<Node> end, std::shared_ptr<Edge> pred)                             : start(start), weight(weight), end(end), pred(pred),    succ(nullptr) {}
    Edge(std::shared_ptr<Node> start, uint32_t weight, std::shared_ptr<Node> end, std::shared_ptr<Edge> pred, std::shared_ptr<Edge> succ) : start(start), weight(weight), end(end), pred(pred),    succ(succ)    {}
};

// Usefull when starting Dijkstra search and we intend to go back to the beggining node
std::shared_ptr<Edge> constructSelfLoopInitEdge(std::shared_ptr<Node> v);

#endif

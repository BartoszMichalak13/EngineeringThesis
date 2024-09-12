/*
  FILE TO BE REVIEWED 

  potential multithreading in the future
*/

#include <cstdint>
#include <optional> //learn how to use it 
#include <list> //array?

struct NeighbourhoodList
{
  //matrix?
  //list of nodes? ||
  //list of edges? ||
  
};

struct Node
{
  uint16_t id; // think how many we may have 
  bool visited; //potentially optional


  NeighbourhoodList* neighbours; // potentially many implemenatitions 
  // one global instance? I think so


  Node(uint16_t id)
  {
    this->id = id;
    visited = false;

  };
  ~Node(){};
};

class Graph {
  public:
    Graph(uint16_t numberOfNodes, uint16_t numberOfEdges)
    {
      this->numberOfNodes = numberOfNodes;
      this->numberOfEdges = numberOfEdges;
      listOfNodes = new Node[numberOfNodes]; //is this correct?
    };
    
    ~Graph()
    {
      delete[] listOfNodes;
    };

    void switchRepresentation();
    void addNode(uint16_t position, uint16_t id)
    {
      listOfNodes[position].id = id;
    }
    //void printData(); // representation? numbers etc

  private:
    uint16_t numberOfNodes;
    uint16_t numberOfEdges;

    Node* listOfNodes;
    
    //list of nodes
    //list of edges?

};


// Mozemy podac graph jako ref a graph tworzymy gdzie indziej
void generateGraph(
    uint16_t numberOfNodes,
    uint16_t numberOfEdges, //optional?
    float denisty)        //potentially more params
{
  Graph graph = Graph(numberOfNodes, numberOfEdges);
  for (uint16_t i = 0; i < numberOfNodes; ++i) 
  {
    graph.listOfNodes[i] = Node(i); // given IDs to nodes is pointless, we may do it at init, here we add edges
  }

  for (uint16_t i = 0; i < numberOfEdges; ++i) 
  {
    /*
    Simple idae: random 2 nodes, connect, repeart till target is met
    */
    uint16_t random1;
    uint16_t random2;
    //add edge 1-2 2-1 bc biderectional 
    // update neighourhood list?
  }
}


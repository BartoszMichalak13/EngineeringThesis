// #include <cstdint>
// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include <string>
// #include <random> 

// struct Node
// {
//   uint16_t id; // think how many we may have 
//   bool visited; //potentially optional

//   Node(uint16_t id)
//   {
//     this->id = id;
//     visited = false;
//   };
//   ~Node(){};

//   operator int() const {
//     return (int)id;
//   }
//   operator uint16_t() const {
//     return id;
//   }
// };

// class Graph {
//   private:
//     uint16_t numberOfNodes;
//     uint16_t numberOfEdges;

//     //adjacency list - table of vectors (vec of vecs if we allow adding new Nodes)
//     std::vector<Node>* adjacencyList;

//   public:
//     Graph(uint16_t numberOfNodes, uint16_t numberOfEdges)
//     {
//       this->numberOfNodes = numberOfNodes;
//       this->numberOfEdges = numberOfEdges;
//       this->adjacencyList = new std::vector<Node>[numberOfNodes];
//     };
//     ~Graph(){};
//     void addEdge(Node node1, Node node2)
//     {
//       // std::cout << "adding edge between " << (uint16_t) node1 << " and " << (uint16_t) node2 << std::endl;
//       this->adjacencyList[(uint16_t)node1].push_back(node2);
//     }
//     //rm edge/vertex
//     void printData()  // mozna ladniej np wszystkie po 3/4 cyfry itd
//     {
//       for (int i = 0; i < numberOfNodes; ++i)
//       {
//         std::cout << "node "<< i << "\n\t"; 
//         for (int j = 0; j < adjacencyList[i].size(); ++j)
//         {
//           std::cout <<  (uint16_t) adjacencyList[i].at(j) << " "; 
//         }
//         std::cout << std::endl;
//       }
//     }
// };

// class Random {
//   private:
//     std::random_device rd;
//     std::mt19937 generator;
//   public:
//     Random() : generator(rd()){};
//     ~Random(){};
//     uint16_t generateRandomNumber(uint16_t max)
//     {
//       return std::uniform_int_distribution<uint16_t> {0, max}(this->generator);
//     }
// };

// // Mozemy podac graph jako ref a graph tworzymy gdzie indziej
// void generateGraph(
//     uint16_t numberOfNodes,
//     uint16_t numberOfEdges, //optional?
//     float denisty)        //potentially more params
// {
//   Graph graph = Graph(numberOfNodes, numberOfEdges);
//   Random randomGen;
//   uint16_t range = numberOfNodes - 1;// std::max({numberOfNodes,numberOfEdges});
//   for (uint16_t i = 0; i < numberOfEdges; ++i) 
//   {
//     /*
//     Simple idae: random 2 nodes, connect, repeart till target is met
//     */ 
//     graph.addEdge( Node(randomGen.generateRandomNumber(range)), Node(randomGen.generateRandomNumber(range)) );
//     //add edge 1-2 2-1 bc biderectional 
//     // update neighourhood list?
//   }
//   graph.printData();
// }

// void writeOutput()
// {

// }

// void parseInput(
//     int argc,
//     char *argv[],
//     uint16_t &numberOfNodes, 
//     uint16_t &numberOfEdges,
//     uint16_t &state)
// {
//   switch(argc) 
//   {
//     case 1: 
//       std::cerr << "No arguments" << std::endl;
//       state = 1;
//       break;
//     default:
//       numberOfNodes = (uint16_t) std::stoi(argv[1]); 
//       numberOfEdges = (uint16_t) std::stoi(argv[2]); 
//       /*
//       ...
//       */
//       state = 0;
//       break;
//   }
// }

// /*
// Usage: numberOfNodes numberOfEdges
// */
// int main(int argc, char *argv[]) 
// {
//   uint16_t numberOfNodes, numberOfEdges, state;

//   //tmp solution
//   float denisty = 0.1;

//   parseInput(argc, argv, numberOfNodes, numberOfEdges, state);
//   if(state) 
//   {
//     std::cerr << "Returning state = " << state << std::endl;
//     return state;
//   }
//   generateGraph(numberOfNodes, numberOfEdges, denisty);
//   writeOutput();
//   return 0;
// }



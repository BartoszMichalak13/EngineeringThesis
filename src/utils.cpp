#include <iostream>
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "algorithmRunner.hpp"

#include <fstream>
#include <sstream>

// Parse function

// Function to parse number of nodes and edges
void parseGraphProperties(
    const std::string& filename,
    uint32_t& numberOfNodes,
    uint32_t& numberOfEdges)
{
  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    if (line.find("SECTION Graph") != std::string::npos) {
      // Parse nodes and edges
      while (std::getline(file, line) && line != "END") {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;

        if (type == "Nodes") {
          lineStream >> numberOfNodes;
        } else if (type == "Edges") {
          lineStream >> numberOfEdges;
        }
      }
      break;
    }
  }
}

// Function to parse edges and terminals
void parseEdgesAndTerminals(
    const std::string& filename,
    std::shared_ptr<Graph>& graph,
    std::vector<uint32_t>& terminals)
{
  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    if (line.find("SECTION Graph") != std::string::npos) {
      // Skip until edge definitions
      while (std::getline(file, line) && line != "END") {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;

        if (type == "E") {
          uint32_t node1, node2, weight;
          lineStream >> node1 >> node2 >> weight;
          graph->addEdge(node1 - 1, node2 - 1, weight); // Adjust for 0-based indexing
        }
      }
    } else if (line.find("SECTION Terminals") != std::string::npos) {
      // Parse terminals
      while (std::getline(file, line) && line != "END") {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;
        if (type == "T") {
          uint32_t terminal;
          lineStream >> terminal;
          terminals.push_back(terminal - 1); // Adjust for 0-based indexing
        }
      }
    }
  }
}

void printUsage(std::string programName) {
  std::cerr << "Usage 0: " << programName << "<mode> <filename>" << std::endl;
  std::cerr << "Usage 1: " << programName << "<mode> <numberOfNodes> <numberOfEdges>" << std::endl;
}

void printModes() {
  std::cerr << "Modes:" <<std::endl;
  std::cerr << "0 - load graph from file" << std::endl;
  std::cerr << "1 - generate graph with given number of nodes and edges" << std::endl;
}

// Function to parse input and update parameters
std::shared_ptr<Graph> parseInput(
    int argc,
    char *argv[],
    uint32_t &numberOfNodes,
    uint32_t &numberOfEdges,
    uint32_t &state,
    bool &printFlag,
    std::vector<uint32_t> &terminals)
{
  std::shared_ptr<Graph> graph;
  switch(argc) 
  {
    case 1: 
      std::cerr << "No arguments" << std::endl;

      printUsage(argv[0]);
      printModes();
      state = 1;
      // return std::shared_ptr<Graph>(new Graph(0,0,0));
      return dummySharedPointerGraph();
      break;
    default:
      uint8_t mode = static_cast<uint8_t>(std::stoi(argv[1]));
      switch (mode)
      {
        case 0:
            parseGraphProperties(argv[2], numberOfNodes, numberOfEdges);

            graph = std::shared_ptr<Graph>(new Graph(numberOfNodes, numberOfEdges, printFlag));

            parseEdgesAndTerminals(argv[2], graph, terminals);
            state = 0;
            return graph;
          break;
        case 1:
            numberOfNodes = static_cast<uint32_t>(std::stoi(argv[2]));
            numberOfEdges = static_cast<uint32_t>(std::stoi(argv[3]));
            state = 0;
            return generateGraph(numberOfNodes, numberOfEdges, printFlag);
          break;

        default:
            std::cerr << "Wrong mode selected" << std::endl;
            printModes();

            state = 1;
            return dummySharedPointerGraph();
          break;
      }

      break;
  }
}

// Function to output graph data or any relevant information
void writeOutput()
{
  // Placeholder for writing output to a file or console
  std::cout << "Output written successfully" << std::endl;
}

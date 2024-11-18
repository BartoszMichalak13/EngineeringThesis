#include <iostream>
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "helpers.hpp"
#include "algorithmRunner.hpp"

#include <fstream>
#include <sstream>
#include <filesystem>

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

//TODO QOL change it so DreyfusWagner is 0th bit
void printAlgorithmsToRunForm() {
  std::cerr << "algorithmsToRun is number which bits correspond to certain algorithms:" << std::endl;
  std::cerr << "0th bit turned to 1 allows to run TakahashiMatsuyama approximation algorith" << std::endl;
  std::cerr << "1st bit turned to 1 allows to run KouMarkowskyBerman approximation algorith" << std::endl;
  std::cerr << "2nd bit turned to 1 allows to run DreyfusWagner exact algorith" << std::endl;

  std::cerr << "Example: algorithmsToRun = 7 => 0...0111 runs all algorithms" << std::endl;
  std::cerr << "         algorithmsToRun = 5 => 0...0101 runs DreyfusWagner and TakahashiMatsuyama" << std::endl;
  std::cerr << "         algorithmsToRun = 2 => 0...0010 runs KouMarkowskyBerman" << std::endl;
}

void printUsage(std::string programName) {
  std::cerr << "Usage 0: " << programName << " <mode> <algorithmsToRun> <filename>" << std::endl;
  std::cerr << "Usage 1: " << programName
            << " <mode> <numberOfNodes> <numberOfEdges> <algorithmsToRun>" << std::endl;
  std::cerr << "Usage 2: " << programName
            << " <mode> <density> <algorithmsToRun> <startingNumberOfNodes> <targetNumberOfNodes> <fileToSaveResults>"
            << std::endl;
  printAlgorithmsToRunForm();
}

void printModes() {
  std::cerr << "Modes:" <<std::endl;
  std::cerr << "0 - load graph from file" << std::endl;
  std::cerr << "1 - generate graph with given number of nodes and edges" << std::endl;
  std::cerr << "2 - generate multiple graphs with given density and range of |V| and saves results to file" << std::endl;
}

// Function to parse input and update parameters
std::shared_ptr<Graph> parseInput(
    int argc,
    char *argv[],
    std::string &fileName,
    uint32_t &numberOfNodes,
    uint32_t &numberOfEdges,
    float    &density,
    uint32_t &algorithmsToRun,
    uint32_t &startingNumberOfNodes,
    uint32_t &targetNumberOfNodes,
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
      state = 1000;
      // return std::shared_ptr<Graph>(new Graph(0,0,0));
      return dummySharedPointerGraph();
      break;
    default:
      uint8_t mode = static_cast<uint8_t>(std::stoi(argv[1]));
      switch (mode)
      {
        case 0:
            if (argc < 5) {
              state = 1000;
              printUsage(argv[0]);
              printModes();
              return dummySharedPointerGraph();
            }
            algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[2]));
            parseGraphProperties(argv[3], numberOfNodes, numberOfEdges);

            graph = std::shared_ptr<Graph>(new Graph(numberOfNodes, numberOfEdges, printFlag));

            parseEdgesAndTerminals(argv[3], graph, terminals);
            fileName = std::filesystem::current_path();
            fileName += argv[4];
            state = mode;
            return graph;
          break;
        case 1:
            if (argc < 5) {
              state = 1000;
              printUsage(argv[0]);
              printModes();
              return dummySharedPointerGraph();
            }
            numberOfNodes = static_cast<uint32_t>(std::stoi(argv[2]));
            numberOfEdges = static_cast<uint32_t>(std::stoi(argv[3]));
            algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[4]));
            startingNumberOfNodes = 0;
            targetNumberOfNodes = 0;
            state = mode;
            return generateGraph(numberOfNodes, numberOfEdges, printFlag);
          break;
        case 2:
          if (argc < 7) {
            state = 1000;
            printUsage(argv[0]);
            printModes();
            return dummySharedPointerGraph();
          }
          numberOfNodes = 0;
          numberOfEdges = 0;
          density = std::stof(argv[2]);
          algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[3]));
          startingNumberOfNodes = static_cast<uint32_t>(std::stoi(argv[4]));
          targetNumberOfNodes = static_cast<uint32_t>(std::stoi(argv[5]));
          // fileName = std::getenv("HOME");
          // fileName = std::filesystem::current_path();

          // // Get the current working directory
          // std::filesystem::path currentDir = std::filesystem::current_path();

          // // Get the parent directory of the current directory
          // std::filesystem::path parentDir = currentDir.parent_path();
          // // fileName = parentDir;

          fileName = std::filesystem::current_path().parent_path();
          fileName += "/results/";
          fileName += argv[6];
          fileName += ".results";

          state = mode;
          return dummySharedPointerGraph(); //generateGraph(numberOfNodes, numberOfEdges, printFlag);
        break;

        default:
            std::cerr << "Wrong mode selected" << std::endl;
            printModes();

            state = 1000;
            return dummySharedPointerGraph();
          break;
      }

      break;
  }
}

// Function to output graph data or any relevant information
void writeOutput(
    const std::string& fileName,
    uint32_t numberOfTerminals,
    uint32_t numberOfNodes,
    uint32_t numberOfEdges,
    uint32_t DreyfusWagnerCost,
    uint32_t DreyfusWagnerTime,
    uint32_t TakahashiMatsuyamaCost,
    uint32_t TakahashiMatsuyamaTime,
    uint32_t KouMarkowskyBermanCost,
    uint32_t KouMarkowskyBermanTime)
{

  // Extract the directory part of the file path
  std::filesystem::path filePath(fileName);
  std::filesystem::path directory = filePath.parent_path();

  // Create directories if they don't exist
  if (!directory.empty() && !std::filesystem::exists(directory)) {
      try {
          std::filesystem::create_directories(directory);
          std::cout << "Directory created: " << directory << std::endl;
      } catch (const std::filesystem::filesystem_error& e) {
          std::cerr << "Error: Unable to create directory " << directory << ": " << e.what() << std::endl;
          return;
      }
  }

  bool fileExists = std::filesystem::exists(fileName);

  // Open file in append mode
  std::ofstream outFile(fileName, std::ios::app);

  // Check if file opened successfully
  if (!outFile.is_open()) {
    std::cerr << "Error: Unable to open file " << fileName << std::endl;
    return;
  }

  // If file doesn't exist, write the header
  if (!fileExists) {
    outFile << "NumberOfTerminals,NumberOfNodes,NumberOfEdges,"
            << "DreyfusWagnerCost,DreyfusWagnerTime,"
            << "TakahashiMatsuyamaCost,TakahashiMatsuyamaTime,"
            << "KouMarkowskyBermanCost,KouMarkowskyBermanTime\n";
  }

  // Write the data as a new row
  outFile << numberOfTerminals << ","
          << numberOfNodes << ","
          << numberOfEdges << ","
          << DreyfusWagnerCost << ","
          << DreyfusWagnerTime << ","
          << TakahashiMatsuyamaCost << ","
          << TakahashiMatsuyamaTime << ","
          << KouMarkowskyBermanCost << ","
          << KouMarkowskyBermanTime << "\n";

  // Close the file
  outFile.close();

  std::cout << "Results written to " << fileName << std::endl;
}

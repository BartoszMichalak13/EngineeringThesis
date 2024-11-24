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

  while (std::getline(file, line))
  {
    if (line.find("SECTION Graph") != std::string::npos)
    {
      // Parse nodes and edges
      while (std::getline(file, line) && line != "END")
      {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;

        if (type == "Nodes")
        {
          lineStream >> numberOfNodes;
        }
        else if (type == "Edges")
        {
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

  while (std::getline(file, line))
  {
    if (line.find("SECTION Graph") != std::string::npos)
    {
      // Skip until edge definitions
      while (std::getline(file, line) && line != "END")
      {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;

        if (type == "E")
        {
          uint32_t node1, node2, weight;
          lineStream >> node1 >> node2 >> weight;
          graph->addEdge(node1 - 1, node2 - 1, weight); // Adjust for 0-based indexing
        }
      }
    }
    else if (line.find("SECTION Terminals") != std::string::npos)
    {
      // Parse terminals
      while (std::getline(file, line) && line != "END")
      {
        std::istringstream lineStream(line);
        std::string type;
        lineStream >> type;
        if (type == "T")
        {
          uint32_t terminal;
          lineStream >> terminal;
          terminals.push_back(terminal - 1); // Adjust for 0-based indexing
        }
      }
    }
  }
}

//TODO QOL change it so DreyfusWagner is 0th bit
void printAlgorithmsToRunForm()
{
  std::cerr << "algorithmsToRun is number which bits correspond to certain algorithms:" << std::endl;
  std::cerr << "0th bit turned to 1 allows to run TakahashiMatsuyama approximation algorith" << std::endl;
  std::cerr << "1st bit turned to 1 allows to run KouMarkowskyBerman approximation algorith" << std::endl;
  std::cerr << "2nd bit turned to 1 allows to run DreyfusWagner exact algorith" << std::endl;

  std::cerr << "Example: algorithmsToRun = 7 => 0...0111 runs all algorithms" << std::endl;
  std::cerr << "         algorithmsToRun = 5 => 0...0101 runs DreyfusWagner and TakahashiMatsuyama" << std::endl;
  std::cerr << "         algorithmsToRun = 2 => 0...0010 runs KouMarkowskyBerman" << std::endl;
}

void printUsage(std::string programName)
{
  std::cerr << "Usage 0: " << programName << " <mode> <algorithmsToRun> <filename> <fileToSaveResults> <printFlag>" << std::endl;
  std::cerr << "Usage 1: " << programName
            << " <mode> <numberOfNodes> <numberOfEdges> <algorithmsToRun> <printFlag>" << std::endl;
  std::cerr << "Usage 2: " << programName
            << " <mode> <density> <algorithmsToRun> <LowerLimitOfNodes> <UpperLimitOfNodes> <fileToSaveResults> <printFlag>"
            << std::endl;
  printAlgorithmsToRunForm();
}

void printModes()
{
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
    uint32_t &LowerLimitOfNodes,
    uint32_t &UpperLimitOfNodes,
    uint32_t &state,
    bool &printFlag,
    std::vector<uint32_t> &terminals)
{
  std::shared_ptr<Graph> graph;
  switch(argc)
  {
    // No arguments Error
    case 1:
      std::cerr << "No arguments" << std::endl;

      printUsage(argv[0]);
      printModes();
      state = 1000;
      return dummySharedPointerGraph();
      break;
    default:
      uint8_t mode = static_cast<uint8_t>(std::stoi(argv[1]));
      switch (mode)
      {
        // Read from file
        case 0:
          // Check if number of arguments is correct
          if (argc < 6)
          {
            state = 1000;
            printUsage(argv[0]);
            printModes();
            return dummySharedPointerGraph();
          }
          algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[2]));
          printFlag = static_cast<bool>(std::stoi(argv[5]));

          parseGraphProperties(argv[3], numberOfNodes, numberOfEdges);

          graph = std::shared_ptr<Graph>(new Graph(numberOfNodes, numberOfEdges, printFlag));

          parseEdgesAndTerminals(argv[3], graph, terminals);
          fileName = std::filesystem::current_path();
          fileName += argv[4];
          state = mode;
          return graph;
          break;
        // Generate single graph
        case 1:
          // Check if number of arguments is correct
          if (argc < 6)
          {
            state = 1000;
            printUsage(argv[0]);
            printModes();
            return dummySharedPointerGraph();
          }
          numberOfNodes = static_cast<uint32_t>(std::stoi(argv[2]));
          numberOfEdges = static_cast<uint32_t>(std::stoi(argv[3]));
          algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[4]));
          printFlag = static_cast<bool>(std::stoi(argv[5]));
          LowerLimitOfNodes = 0;
          UpperLimitOfNodes = 0;
          fileName = "";
          state = mode;
          return generateGraph(numberOfNodes, numberOfEdges, printFlag);
          break;
        // Generate graphs in loop
        case 2:
          // Check if number of arguments is correct
          if (argc < 8)
          {
            state = 1000;
            printUsage(argv[0]);
            printModes();
            return dummySharedPointerGraph();
          }
          numberOfNodes = 0;
          numberOfEdges = 0;
          density = std::stof(argv[2]);
          algorithmsToRun = static_cast<uint32_t>(std::stoi(argv[3]));
          LowerLimitOfNodes = static_cast<uint32_t>(std::stoi(argv[4]));
          UpperLimitOfNodes = static_cast<uint32_t>(std::stoi(argv[5]));

          fileName = std::filesystem::current_path().parent_path();
          fileName += "/results/";
          fileName += argv[6];
          fileName += ".results";
          printFlag = static_cast<bool>(std::stoi(argv[7]));

          state = mode;
          return dummySharedPointerGraph(); //generateGraph(numberOfNodes, numberOfEdges, printFlag);
          break;

        // Wrong mode selected
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

void writeOutput(
    const std::string& fileName,
    uint32_t numberOfTerminals,
    uint32_t numberOfNodes,
    uint32_t numberOfEdges,
    uint32_t DreyfusWagnerCost,
    uint32_t DreyfusWagnerTime,
    uint32_t TakahashiMatsuyamaCost,
    uint32_t TakahashiMatsuyamaTime,
    std::vector<std::chrono::microseconds> TMtimes,
    uint32_t KouMarkowskyBermanCost,
    uint32_t KouMarkowskyBermanTime,
    std::vector<std::chrono::microseconds> KMBtimes
    )
{
  // Extract the directory part of the file path
  std::filesystem::path filePath(fileName);
  std::filesystem::path directory = filePath.parent_path();

  // Create directories if they don't exist
  if (!directory.empty() && !std::filesystem::exists(directory))
  {
    try
    {
      std::filesystem::create_directories(directory);
      std::cout << "Directory created: " << directory << std::endl;
    }
    catch (const std::filesystem::filesystem_error& e)
    {
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
            << "TMdurationNotNecessary,TMdurationInit,TMdurationPrepareForNextIteration,TMdurationMainLoop,"
            << "KouMarkowskyBermanCost,KouMarkowskyBermanTime,"
            << "KMBdurationNotNecessary,KMBdurationStep1,KMBdurationStep2,"
            << "KMBdurationStep3,KMBdurationStep4,KMBdurationStep5\n";
  }

  // Write the data as a new row
  outFile << numberOfTerminals << ","
          << numberOfNodes << ","
          << numberOfEdges << ","
          << DreyfusWagnerCost << ","
          << DreyfusWagnerTime << ","
          << TakahashiMatsuyamaCost << ","
          << TakahashiMatsuyamaTime << ","
          << TMtimes.at(0).count() << ","     // TMdurationNotNecessary
          << TMtimes.at(1).count() << ","     // TMdurationInit
          << TMtimes.at(2).count() << ","     // TMdurationPrepareForNextIteration
          << TMtimes.at(3).count() << ","     // TMdurationMainLoop
          << KouMarkowskyBermanCost << ","
          << KouMarkowskyBermanTime << ","
          << KMBtimes.at(0).count() << ","    // KMBdurationNotNecessary
          << KMBtimes.at(1).count() << ","    // KMBdurationStep1
          << KMBtimes.at(2).count() << ","    // KMBdurationStep2
          << KMBtimes.at(3).count() << ","    // KMBdurationStep3
          << KMBtimes.at(4).count() << ","    // KMBdurationStep4
          << KMBtimes.at(5).count() << "\n";  // KMBdurationStep5


  // Close the file
  outFile.close();

  std::cout << "Results written to " << fileName << std::endl;
}

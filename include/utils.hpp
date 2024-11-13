#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>

// Function declarations
// std::shared_ptr<Graph> parseInput(
//     int argc,
//     char *argv[],
//     uint32_t &numberOfNodes,
//     uint32_t &numberOfEdges,
//     uint32_t &state,
//     bool &printFlag,
//     std::vector<uint32_t> terminals);

// Function to parse input and update parameters
std::shared_ptr<Graph> parseInput(
    int argc,
    char *argv[],
    uint32_t &numberOfNodes,
    uint32_t &numberOfEdges,
    uint32_t &state,
    bool &printFlag,
    std::vector<uint32_t> &terminals);

void writeOutput();


// #include <iostream>
// #include "graph.hpp"
// #include "random.hpp"
// #include "utils.hpp"
// #include "algorithmRunner.hpp"

// #include <fstream>
// #include <sstream>

// Parse function

// // Function to parse number of nodes and edges
// void parseGraphProperties(
//     const std::string& filename,
//     uint32_t& numberOfNodes,
//     uint32_t& numberOfEdges);

// // Function to parse edges and terminals
// void parseEdgesAndTerminals(
//     const std::string& filename,
//     std::shared_ptr<Graph>& graph,
//     std::vector<uint32_t>& terminals)

// void printUsage(std::string programName) {
//   std::cerr << "Usage 0: " << programName << "<mode> <filename>" << std::endl;
//   std::cerr << "Usage 1: " << programName << "<mode> <numberOfNodes> <numberOfEdges>" << std::endl;
// }

// void printModes() {
//   std::cerr << "Modes:" <<std::endl;
//   std::cerr << "0 - load graph from file" << std::endl;
//   std::cerr << "1 - generate graph with given number of nodes and edges" << std::endl;
// }


#endif

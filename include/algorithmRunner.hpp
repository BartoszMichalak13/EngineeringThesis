#ifndef ALGORITHM_RUNNER_HPP
#define ALGORITHM_RUNNER_HPP

#include <cstdint>
#include <chrono>

std::shared_ptr<Graph> generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool  printFlag);
void runAlgorithms(std::shared_ptr<Graph> graph, std::vector<uint32_t> terminals, uint32_t algorithmsToRun, const std::string& fileName);

#endif

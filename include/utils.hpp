#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>

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
    std::vector<uint32_t> &terminals);

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
    uint32_t KouMarkowskyBermanTime);


#endif

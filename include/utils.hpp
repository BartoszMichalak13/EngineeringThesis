#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>

// TODO ta nazwa kłamie - split into multiple functions?
/*
TODO fill description
*/
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

/*
Writes all passed arguments to a file. If it or it's directories doesn't exist it should create it.
*/
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
    std::vector<std::chrono::microseconds> KMBtimes);

#endif

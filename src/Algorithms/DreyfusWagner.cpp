#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <queue>
#include <vector>
#include <memory>
#include <limits>
#include <iomanip>

#include "graph.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"

bool DreyfusPrintFlag = false;

bool compareUintMatrices(std::vector<std::vector<uint32_t>> matrix1, std::vector<std::vector<uint32_t>> matrix2, uint32_t size) {
  bool status = true;
  // Check if matrices have the specified size
  if (matrix1.size() != size || matrix2.size() != size) {
    std::cerr << "Matrices do not match the specified size." << std::endl;
    status = false;
  }

  for (uint32_t i = 0; i < size; ++i) {
    // Check if each row has the correct number of columns
    if (matrix1[i].size() != size || matrix2[i].size() != size) {
      std::cerr << "Row " << i << " does not match the specified size." << std::endl;
      status = false;
    }

    for (uint32_t j = 0; j < size; ++j) {
      if (matrix1[i][j] != matrix2[i][j]) {
        std::cout << "Difference at position (" << i << ", " << j << "): "
                  << "matrix1 = " << matrix1[i][j] << ", "
                  << "matrix2 = " << matrix2[i][j] << std::endl;
        status = false;
      }
    }
  }
  if (status)
    std::cout << "Matrices are identical." << std::endl;

  return status;
}

void printVectorOfSets(std::vector<std::set<uint32_t>>& vecotrOfSets){
  for (auto const &set : vecotrOfSets) {
    std::cout << "Subset = ";
    for (auto setItem : set)
      std::cout << setItem << " ";
    std::cout << "\n";
  }
}

bool findInSet(std::set<uint32_t> A, uint32_t value){
  for (uint32_t setItem : A)
    if (setItem == value)
      return true;
  return false;
}

std::map<std::set <uint32_t>, uint32_t> subsetMap;

std::vector<std::set<uint32_t>> generateSubsets(std::vector<uint32_t> terminals) {
  const uint32_t terminalsSize = terminals.size();
  const uint32_t count = pow(2, terminalsSize);

  std::vector<std::set<uint32_t>> allSubsets;
  uint32_t subsetIndex = 0;

  for (uint32_t i = 0; i < count; ++i)
  {
    std::set<uint32_t> newSubset;
    for (uint32_t j = 0; j < terminalsSize; ++j)
    {
      if ((i & (1 << j)) > 0)
      {
        newSubset.insert(terminals[j]);
      }
    }
    if (newSubset.size() > 1)
    {
      allSubsets.push_back(newSubset);
      subsetMap.insert(make_pair(newSubset,subsetIndex));
      subsetIndex++;
    }
  }

  return allSubsets;
}

std::set<uint32_t> getSetFromMap(uint32_t idx)
{
  std::set<uint32_t> subset;
  for(std::map<std::set<uint32_t>,uint32_t >::const_iterator it = subsetMap.begin(); it != subsetMap.end(); ++it)
  {
    if (it -> second == idx)
    {
      subset = it->first;
    }
  }
  return subset;
}

void initializeSteinerTable(
    std::vector<std::vector<uint32_t>> &steinerDistanceTable,
    std::vector<std::vector<uint32_t>> &graph,
    uint32_t numberOfNodes,
    uint32_t sizeOfTable)
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
  {
    for (uint32_t j = 0; j < numberOfNodes; ++j)
    {
      steinerDistanceTable[i][j] = graph[i][j];
    }
  }
  for (uint32_t i = numberOfNodes; i < sizeOfTable; ++i) {
    for (uint32_t j = numberOfNodes; j < sizeOfTable; ++j) {
      steinerDistanceTable[i][j] = std::numeric_limits<uint32_t>::max();
      steinerDistanceTable[j][i] = std::numeric_limits<uint32_t>::max();
      if (j == i) {
        steinerDistanceTable[i][j] = 0;
      }
    }
  }

}

uint32_t fetchIndexofMapofSets(std::set<uint32_t> subset)
{
  bool wasInit = false;
  uint32_t key;
  for(std::map<std::set<uint32_t>,uint32_t >::const_iterator it = subsetMap.begin(); it != subsetMap.end(); ++it)
  {
    if (it -> first == subset)
    {
      key = it ->second;
      wasInit = true;
    }
  }
  if (wasInit)
    return key;
  else {
    std::cout << "No such a subset found in fetchIndexofMapofSets \n";
    std::cout << "{ ";
    for (const auto& elem : subset) {
        std::cout << elem << " ";
    }
    std::cout << "}" << std::endl;
    return 0;
  }
}

/*
  If you have a std::vector<std::vector<int>> G and you call G.size(),
  it will return the number of rows (or outer vectors) in G - by GPT
*/
uint32_t Graph::calculateSteiner(
    std::vector<uint32_t> C,
    std::vector<std::vector<uint32_t>> adjMatrix,
    std::vector<std::set<uint32_t>> allSubsets,
    uint32_t q)
{
  uint32_t sizetable = allSubsets.size() + numberOfNodes;
  std::vector<std::vector<uint32_t>> steinerDistanceTable(
      sizetable,
      std::vector<uint32_t> (sizetable, std::numeric_limits<uint32_t>::max()));

  //initialise steiner table
  initializeSteinerTable(steinerDistanceTable, adjMatrix, numberOfNodes, sizetable);

  uint32_t indexOfD = 0;

  std::set<std::set<uint32_t>>::iterator setIterator; // iterator for the "outer" structure

  constexpr uint32_t ONE = 1;

  for (uint32_t m = 2; m <= C.size(); ++m)
  {
    // for each subset
    for (uint32_t i = 0; i < allSubsets.size(); ++i)
    {
      std::set<uint32_t> D = allSubsets[i];
      if (D.size() == m)
      {
        indexOfD = subsetMap[D] + numberOfNodes;

        // initialise distance from elementOfD to every node in graph
        for (uint32_t k = 0; k < numberOfNodes; ++k) {
          steinerDistanceTable[indexOfD][k] = std::numeric_limits<uint32_t>::max();
          steinerDistanceTable[k][indexOfD] = std::numeric_limits<uint32_t>::max();
        }

        // for every element in graph
        for (uint32_t j = 0; j < numberOfNodes; ++j)
        {
          uint32_t u = std::numeric_limits<uint32_t>::max(); ;

          std::vector<uint32_t> elements(D.begin(), D.end());
          uint32_t n = elements.size();

          // Iterate over all possible subsets (2^n subsets)
          for (uint32_t mask = 0; mask < (ONE << n); ++mask) {
            std::set<uint32_t> E;
            for (uint32_t k = 0; k < n; ++k) {
              if (mask & (ONE << k)) { // Check if the i-th bit is set
                E.insert(elements[k]);
              }
            }
            uint32_t indexOfE;// = fetchIndexofMapofSets(E) + numberOfNodes;
            if(E.size() == 0) //empty set
              continue;
            if(E.size() == 1)
            {
              std::set<uint32_t>::iterator node = E.begin();
              indexOfE = *node;
            }
            else
            {
              indexOfE = fetchIndexofMapofSets(E) + numberOfNodes;
            }
            uint32_t distanceFromEToJ = steinerDistanceTable[indexOfE][j];

            std::set<uint32_t> DMinusE;
            DMinusE = D;

            for (std::set<uint32_t>::iterator iter = DMinusE.begin(); iter != DMinusE.end();)
            {
              if (E.find(*iter) != E.end())
              {
                iter = DMinusE.erase(iter);
              }
              else
              {
                ++iter;
              }
            }
            uint32_t dist, indexOfSubsetDMinusE;
            if(DMinusE.size() == 0)
              continue;
            if(DMinusE.size() == 1)
            {
              std::set<uint32_t>::iterator node = DMinusE.begin();
              indexOfSubsetDMinusE = *node;
            }
            else
            {
              indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE) + numberOfNodes;
            }

            uint32_t distanceFromDMinusEToJ = steinerDistanceTable[indexOfSubsetDMinusE][j];
            if (distanceFromEToJ < std::numeric_limits<uint32_t>::max() && distanceFromDMinusEToJ < std::numeric_limits<uint32_t>::max()) {
              dist = distanceFromEToJ + distanceFromDMinusEToJ;
            } else {
              dist = std::numeric_limits<uint32_t>::max();
            }
            u = std::min(u, dist);
          }

          uint32_t val;

          for (uint32_t k = 0; k < numberOfNodes; ++k)
          {
            val = adjMatrix[k][j];
            if (val < std::numeric_limits<uint32_t>::max() && u < std::numeric_limits<uint32_t>::max()) {
              if (val + u < steinerDistanceTable[indexOfD][k]) {
                steinerDistanceTable[indexOfD][k] = val + u;
                steinerDistanceTable[k][indexOfD] = val + u;
              }
            }
          }
        }
      }
    }
    if (DreyfusPrintFlag) {
      std::cout << "M's " << m << "/" << C.size() << std::endl;
    }

  }

  uint32_t v = std::numeric_limits<uint32_t>::max();

  for (uint32_t j = 0; j < numberOfNodes; j++)
  {
    uint32_t u = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> elements(C.begin(), C.end());
    uint32_t n = elements.size();

    // Iterate over all possible subsets (2^n subsets)
    for (uint32_t mask = 0; mask < (ONE << n); ++mask) {
      std::set<uint32_t> E;
      for (uint32_t k = 0; k < n; ++k) {
        if (mask & (ONE << k)) { // Check if the i-th bit is set
          E.insert(elements[k]);
        }
      }

      uint32_t indexOfE;
      if(E.size() == 0) //empty set
        continue;
      if(E.size() == 1)
      {
        std::set<uint32_t>::iterator node = E.begin();
        indexOfE = *node;
      }
      else
      {
        indexOfE = fetchIndexofMapofSets(E) + numberOfNodes;
      }
      uint32_t distEJ = steinerDistanceTable[indexOfE][j];

      std::set<uint32_t> subsetCMinusE (C.begin(), C.end());
      for (std::set<uint32_t>::iterator iter = subsetCMinusE.begin(); iter != subsetCMinusE.end();)
      {
        if (E.find(*iter) != E.end())
        {
          iter = subsetCMinusE.erase(iter);
        }
        else
        {
          ++iter;
        }
      }
      uint32_t indexOfSubsetCMinusE;
      if(subsetCMinusE.size() == 0)
        continue;
      if(subsetCMinusE.size() == 1)
      {
        std::set<uint32_t>::iterator node = subsetCMinusE.begin();
        indexOfSubsetCMinusE = *node;
      }
      else
      {
        indexOfSubsetCMinusE = fetchIndexofMapofSets(subsetCMinusE) + numberOfNodes;
      }

      uint32_t val = steinerDistanceTable[indexOfSubsetCMinusE][j];
      if (distEJ < std::numeric_limits<uint32_t>::max() && val < std::numeric_limits<uint32_t>::max()) {

        if (u > distEJ + val) {
          u = distEJ + val;
        }
      }
    }
    if (DreyfusPrintFlag) {
      std::cout << "Final J's " << j << "/" << numberOfNodes << std::endl;
    }

    uint32_t value = adjMatrix[q][j];
    if (allSubsets[allSubsets.size() - 1].find(j) != allSubsets[allSubsets.size() - 1].end()) // j is one of the terminals
    {
      value = 0;
    }


    if (value < std::numeric_limits<uint32_t>::max() && u < std::numeric_limits<uint32_t>::max())
    {
      v = std::min(v, value + u);
    }
    v = std::min(v, steinerDistanceTable[sizetable - 1][j]); //mimal value for any j and all terminals e.i. optimal steiner tree

  }
  if (printFlag) {
    std::cout << " //// ";
    for(uint32_t j = 0; j < sizetable ; j++)
      std::cout << " " << std::setw(4) << std::setfill('0')  << j;
    std::cout << std::endl;

    for(uint32_t i = 0 ; i < sizetable ; i++ )
    {
      std::cout << " " << std::setw(4) << std::setfill('0') << i << " ";
      for(uint32_t j = 0; j < sizetable ; j++)
      {
        if (steinerDistanceTable[i][j] == std::numeric_limits<uint32_t>::max())
        {
          std::cout << " MAXX";
        } else {
          std::cout << " " << std::setw(4) << std::setfill('0')  << steinerDistanceTable[i][j];
        }
      }
      std::cout << " " << std::endl;
    }
    std::cout << "Minimum Steiner Distance " << v << "\n";
  }
  return v;
}

uint32_t Graph::DreyfusWagner(std::vector<uint32_t> terminals) {
  // get graph instance
  std::shared_ptr<Graph> self = shared_from_this();

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  uint32_t q = terminals[0];

  // terminals.erase(terminals.begin());

  std::vector<std::set<uint32_t>> allSubsets = generateSubsets(terminals);
  // terminals.erase(terminals.begin());


  // convert vertices to uint vec
  std::shared_ptr<std::vector<uint32_t>> uintVertices = nodeArrayToUint(vertices, numberOfNodes);

  // calculate all shortest paths
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> tmpShortestPaths = AllPairsShortestPath((*uintVertices));

  std::vector<std::vector<uint32_t>> adjMatrix2 = toAdjacencyMatrix();

  //fill max fields with shortest paths
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    for (uint32_t j = 0; j < numberOfNodes; ++j) {
      uint32_t pathWeight = findShortestPathAndReturnWeight(tmpShortestPaths, i, j);
      if (pathWeight == std::numeric_limits<uint32_t>::max()){
        std::cerr << "DW did not found path from " << i << " to " << j << std::endl;
        return std::numeric_limits<uint32_t>::max();
      }

      if (i != j) {
        if (adjMatrix2[i][j] >= pathWeight) {
          adjMatrix2[i][j] = std::min(adjMatrix2[i][j], pathWeight);
          adjMatrix2[j][i] = std::min(adjMatrix2[j][i], pathWeight);
        } else {
          std::cout << "Path from " << i << " to " << j << " = " << pathWeight << " is longer than " << adjMatrix2[i][j] << std::endl;
        }
      }
    }
  }

  uint32_t steinerWeight = calculateSteiner(
    terminals,
    adjMatrix2,
    allSubsets,
    q);
  subsetMap.clear();

  return steinerWeight;
}
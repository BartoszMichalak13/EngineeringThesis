// https://bitbucket.org/NachiketDan/dreyfus-wagner/src/master/

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <queue>
#include <vector>
#include <memory>
#include <limits>


#include "graph.hpp"
#include "helpers.hpp"
#include "debugPrints.hpp"


#include <iomanip>

//TODO hipoteza czemu wyniki zle: wyniki są zawyżone przez fakt, że biorąc shortest path miedzy 2 wierzcholkami nie zapisujemy krawedzi,
// ergo one moga byc liczne 2x
//
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



// Recursive function to generate subsets
// void generateProperSubsets(const std::set<uint32_t>& fullSet,
//                            std::set<uint32_t>& currentSubset,
//                            std::vector<std::set<uint32_t>>& allProperSubsets,
//                            std::set<uint32_t>::iterator iter) {
//     if (iter == fullSet.end()) {
//         if (!currentSubset.empty() && currentSubset != fullSet) {
//             allProperSubsets.push_back(currentSubset);
//         }
//         return;
//     }

//     // Exclude current element
//     generateProperSubsets(fullSet, currentSubset, allProperSubsets, std::next(iter));

//     // Include current element
//     currentSubset.insert(*iter);
//     generateProperSubsets(fullSet, currentSubset, allProperSubsets, std::next(iter));
//     currentSubset.erase(*iter);
// }

// std::vector<std::set<uint32_t>> getAllProperSubsets(const std::set<uint32_t>& fullSet) {
//     std::vector<std::set<uint32_t>> allProperSubsets;
//     std::set<uint32_t> currentSubset;
//     generateProperSubsets(fullSet, currentSubset, allProperSubsets, fullSet.begin());
//     return allProperSubsets;
// }

//TODO clear this file
//TODO check if floyd is maybe faster?
//TODO how about aadjacency matrix instead of list,
//maybe I'll get extra time for no significant increase in memory usage



std::map<std::set <uint32_t>, uint32_t> subsetMap;



std::vector<std::set<uint32_t>> generateSubsets(std::vector<uint32_t> terminals) { //, uint32_t terminalsSize) { //, uint32_t graphsize) {
  const uint32_t terminalsSize = terminals.size();
  const uint32_t count = pow(2, terminalsSize);

  std::vector<std::set<uint32_t>> allSubsets;
  uint32_t subsetIndex = 0;
  // The outer for loop will run 2^n times to print all subset .
  // Here variable i will act as a binary counter

  for (uint32_t i = 0; i < count; ++i)
  {
    // The inner for loop will run n times , As the maximum number of elements a set can have is n
    // This loop will generate a subset
    std::set<uint32_t> newSubset;
    for (uint32_t j = 0; j < terminalsSize; ++j)
    {
      // This if condition will check if jth bit in binary representation of  i  is set or not
      // if the value of (i & (1 << j)) is greater than 0 , include arr[j] in the current subset
      // otherwise exclude arr[j]
      if ((i & (1 << j)) > 0)
      {
          newSubset.insert(terminals[j]);
      }
    }

    if (newSubset.size() > 1)
    {
      allSubsets.push_back(newSubset);
      //subsetMap[newSubset] = subsetIndex;
      subsetMap.insert(make_pair(newSubset,subsetIndex));
      std::cout << "Adding subset with index" << subsetIndex << std::endl;
      subsetIndex++;
    }
  }

  for (auto const &x : subsetMap)
  {
    std::cout << "Subset = ";
    std::set<uint32_t> element = x.first;
    for (auto f : element)
    {
      std::cout << f << " ";
    }

    std::cout << std::endl;
    std::cout << "Index = " << x.second << std::endl;
  }

  std::cout<<"Completed subsets"<< std::endl;

  return allSubsets;
}










/*
repair
*/
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
  uint32_t key;
  for(std::map<std::set<uint32_t>,uint32_t >::const_iterator it = subsetMap.begin(); it != subsetMap.end(); ++it)
  {
    if (it -> first == subset)
    {
      key = it ->second;
    }
  }

  return key;
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
  // std::shared_ptr<Graph> self = shared_from_this();

  uint32_t sizetable = allSubsets.size() + numberOfNodes; //+ graph.size(); //
  std::vector<std::vector<uint32_t>> steinerDistanceTable(
      sizetable,
      std::vector<uint32_t> (sizetable, std::numeric_limits<uint32_t>::max()));

  std::cout << "works 2" << std::endl;

  //initialise steiner table
  initializeSteinerTable(steinerDistanceTable, adjMatrix, numberOfNodes, sizetable);
  std::cout << "works 3" << std::endl;

  // printMatrix(steinerDistanceTable);
  uint32_t indexOfSubset = 0;

  std::set<std::set<uint32_t>>::iterator setIterator; // iterator for the "outer" structure


  //is it right range?
  for (uint32_t m = 2; m <= C.size(); ++m)
  {


    // for each subset
    for (uint32_t i = 0; i < allSubsets.size(); ++i)
    {
      std::set<uint32_t> elementOfD = allSubsets[i];
      if (elementOfD.size() == m)
      {
        indexOfSubset = subsetMap[elementOfD] + numberOfNodes;
        std::cout << " indexOfSubset: " << indexOfSubset << std::endl;

        // initialise distance from elementOfD to every node in graph
        for (uint32_t k = 0; k < numberOfNodes; ++k) {
          steinerDistanceTable[indexOfSubset][k] = std::numeric_limits<uint32_t>::max();
          steinerDistanceTable[k][indexOfSubset] = std::numeric_limits<uint32_t>::max();
        }

        // for every element in graph
        for (uint32_t j = 0; j < numberOfNodes; ++j)
        {
          uint32_t u = std::numeric_limits<uint32_t>::max(); ;


          //maybe I'm not right but iterates over emenets of size 1
          // std::set<uint32_t>::iterator nodeIterator;
          for (std::set<uint32_t>::iterator nodeIterator = elementOfD.begin(); nodeIterator != elementOfD.end(); ++nodeIterator)
          {
            // E is being considered as a single number from the set named elementOfD
            uint32_t E = *nodeIterator;
            uint32_t distanceFromEToJ = steinerDistanceTable[E][j];

            std::set<uint32_t> DMinusE;
            DMinusE = elementOfD;

            for (std::set<uint32_t>::iterator iter = DMinusE.begin(); iter != DMinusE.end();)
            {
              if (*iter == E)
              {
                iter = DMinusE.erase(iter);
              }
              else
              {
                ++iter;
              }
            }
            uint32_t dist, indexOfSubsetDMinusE;
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
            // std::cout<<"distanceFromDMinusEToJ = "<< distanceFromDMinusEToJ <<"\n";
            // std::cout<<"distanceFromEToJ = "<< distanceFromEToJ <<"\n";
            if (distanceFromEToJ < std::numeric_limits<uint32_t>::max() && distanceFromDMinusEToJ < std::numeric_limits<uint32_t>::max()) {
              dist = distanceFromEToJ + distanceFromDMinusEToJ;
            } else {
              dist = std::numeric_limits<uint32_t>::max();
            }
            u = std::min(u, dist);
          }

          uint32_t indexOfD = fetchIndexofMapofSets(elementOfD) + numberOfNodes;
          std::cout << "indexOfC  "<<indexOfD << std::endl;
          uint32_t value = steinerDistanceTable[indexOfD][j];
          u = std::min(u, value);

          uint32_t val;

          for (uint32_t k = 0; k < numberOfNodes; ++k)
          {
            val = adjMatrix[k][j];



            // if (val != 0)
            // {


              if (val < std::numeric_limits<uint32_t>::max() && u < std::numeric_limits<uint32_t>::max()) {
                if (val + u < steinerDistanceTable[indexOfSubset][k]) {
                    steinerDistanceTable[indexOfSubset][k] = val + u;
                    steinerDistanceTable[k][indexOfSubset] = val + u;
                //     // prevNode[indexOfSubset][k] = j/* Node contributing to the optimal distance */;
                //     prevNode.push_back(j);
                }
                // steinerDistanceTable[indexOfSubset][k] = std::min(steinerDistanceTable[indexOfSubset][k], val + u);
                // steinerDistanceTable[k][indexOfSubset] = std::min(steinerDistanceTable[k][indexOfSubset], val + u);



              // }

              // std::cout << "I is" << k << std::endl;
              // std::cout << "Steiner table data is is" << steinerDistanceTable[indexOfSubset][k] << std::endl;
            }
          }
        }
      }
    }
  }
  // uint32_t vm = std::numeric_limits<uint32_t>::max();
  uint32_t v = std::numeric_limits<uint32_t>::max();

  for (uint32_t j = 0; j < numberOfNodes; j++)
  {
    uint32_t u = std::numeric_limits<uint32_t>::max();

    std::set <uint32_t> C1 (C.begin(),C.end());

      std::cout << "C1 Subset = ";
      for (auto setItem : C1)
        std::cout << setItem << " ";
      std::cout << "\n";
    for (std::set<uint32_t>::iterator setiter = C1.begin(); setiter != C1.end(); ++setiter)
    {
      uint32_t E = *setiter;
      // uint32_t distEJ = adjMatrix[E][j];
      uint32_t distEJ = steinerDistanceTable[E][j];


      std::set<uint32_t> subsetCMinusE (C.begin(), C.end());
      // uint32_t distEJ = adjMatrix[E][j];
      for (std::set<uint32_t>::iterator iter = subsetCMinusE.begin(); iter != subsetCMinusE.end();)
      {
        if (*iter == E)
        {
          iter = subsetCMinusE.erase(iter);
        }
        else
        {
          ++iter;
        }
      }
      uint32_t indexOfSubsetCMinusE;
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
      // std::cout << std::endl;

      // std::cout << " Val bef u: " << val << std::endl;
      // std::cout << "E is: " << E;
      // std::cout << "   J: "<< j << std::endl;
      // std::cout << " distEJ bef u: " << distEJ << std::endl;
      // std::cout << " IndexSubset: " << indexOfSubsetCMinusE;


      //TODO DISTANCE ZBIORU 1 2 3 DO 1 2 CZY 3 > 0?

      // std::cout << std::endl;

      if (distEJ < std::numeric_limits<uint32_t>::max() && val < std::numeric_limits<uint32_t>::max()) {
        std::cout << "CURR BEST U     distEJ (E= "<<E<<" ,J= "<<j<<") is: " << distEJ << std::endl;
        std::cout << "CURR BEST U     val (indexOfSubsetCMinusE= "<<indexOfSubsetCMinusE<<" ,J= "<<j<<")is: " << val << std::endl;
        u = std::min(u, distEJ + val);
        std::cout << "CURR BEST U     u "<<u<< std::endl;

      }

      // std::cout << " U is: " <<u << std::endl;
    }
    //TODO TEGO W ALG NIE BYLO
    uint32_t indexOfC = fetchIndexofMapofSets(C1) + numberOfNodes;
    std::cout << "indexOfC  "<<indexOfC << std::endl;
    uint32_t val = steinerDistanceTable[indexOfC][j];
    u = std::min(u, val);

    // std::cout << " J is: " << j << std::endl;
    // // std::cout << " U final is:" <<u;
    uint32_t value = adjMatrix[q][j];
    // std::cout << "value (D(q,J)) is: " << value << std::endl;
    // std::cout << "value u is: " << u << std::endl;
    // std::cout << "V is: " << v << std::endl;
    // if (value < std::numeric_limits<uint32_t>::max() && value > 0 && u < std::numeric_limits<uint32_t>::max())
    if (value < std::numeric_limits<uint32_t>::max() && u < std::numeric_limits<uint32_t>::max())
    {
      std::cout << "CURR BEST V       (D("<<q<<","<<j<<")) is: " << value << std::endl;
      std::cout << "CURR BEST V       u is: " << u << std::endl;
      // v = std::min(v, value + u);






      v = std::min(v, value + u);






    }
    std::cout << " V is: " << v << std::endl;

  }

  std::cout << " //// ";
  for(uint32_t j = 0; j < sizetable ; j++)
    std::cout << " " << std::setw(4) << std::setfill('0')  << j;
  std::cout << std::endl;

  // for(uint32_t i = 0 ; i < steinerDistanceTable.size() ; i++ )
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
    // std::cout << std::numeric_limits<uint32_t>::max() << std::endl;
    // std::cout << std::numeric_limits<uint32_t>::infinity() << std::endl;
  }

  // printMatrix(adjMatrix);

  // printNodeVector(prevNode);

  std::cout << "Minimum Steiner Distance " << v << "\n";
  return v;
}








void FloydAllPairsShortestPath(std::vector<std::vector<uint32_t>>& a, uint32_t n)
{
  uint32_t k, i, j;
  for (k = 0; k < n; k++)
  {
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        if (a.at(i).at(j) > a.at(i).at(k) + a.at(k).at(j) )
        {
          a.at(i).at(j) = a.at(i).at(k) + a.at(k).at(j);
          // p[i][j] = k;
        }
      }
    }
  }
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
  // std::vector<std::set<uint32_t>> allSubsets = generateSubsets(*uintVertices);

  std::cout << "A -1 " << std::endl;

  printNodeVector(terminals);
  terminals.erase(terminals.begin());
  printNodeVector(terminals);


  // convert vertices to uint vec
  std::shared_ptr<std::vector<uint32_t>> uintVertices = verticesToUint(vertices, numberOfNodes);
  std::cout << "A 0 " << std::endl;


  //TODO czasem tu staje
  // calculate all shortest paths
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> tmpShortestPaths = AllPairsShortestPath((*uintVertices));

  std::cout << "A 1 " << std::endl;

  std::vector<std::vector<uint32_t>> adjMatrix1 = toAdjacencyMatrix();
  FloydAllPairsShortestPath(adjMatrix1, numberOfNodes);
  // std::cout << "A 2 " << std::endl;

  // //fill max fields with shortest paths
  // for (uint32_t i = 0; i < numberOfNodes; ++i) {
  //   for (uint32_t j = 0; j < numberOfNodes; ++j) {
  //     if (adjMatrix1[i][j] == std::numeric_limits<uint32_t>::max()) {
  //       uint32_t pathWeight = findShortestPathAndReturnWeight(tmpShortestPaths, i, j);
  //       if (pathWeight == std::numeric_limits<uint32_t>::max()){
  //         std::cerr << "DW did not found path from " << i << " to " << j << std::endl;
  //         return std::numeric_limits<uint32_t>::max();
  //       }
  //     //  std::cout << "Path weight from " << i << " to " << j << " = " << pathWeight << std::endl;
  //       adjMatrix1[i][j] = pathWeight;
  //       adjMatrix1[j][i] = pathWeight;
  //     }
  //   }
  // }

  std::vector<std::vector<uint32_t>> adjMatrix2 = toAdjacencyMatrix();
  std::cout << "A 2 " << std::endl;

  //fill max fields with shortest paths
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    for (uint32_t j = 0; j < numberOfNodes; ++j) {
      // if (adjMatrix2[i][j] == std::numeric_limits<uint32_t>::max()) {
        uint32_t pathWeight = findShortestPathAndReturnWeight(tmpShortestPaths, i, j);
        if (pathWeight == std::numeric_limits<uint32_t>::max()){
          std::cerr << "DW did not found path from " << i << " to " << j << std::endl;
          return std::numeric_limits<uint32_t>::max();
        }
      //  std::cout << "Path weight from " << i << " to " << j << " = " << pathWeight << std::endl;
      //TODO sometimes minal path is longer than one edge? HOW?
      if (i != j) {
        // std::cout << "Bef from " << i << " to " << j << " = " << adjMatrix2[i][j] << std::endl;
        if (adjMatrix2[i][j] >= pathWeight) {
          adjMatrix2[i][j] = std::min(adjMatrix2[i][j], pathWeight);
          adjMatrix2[j][i] = std::min(adjMatrix2[j][i], pathWeight);
        } else {
          std::cout << "Path from " << i << " to " << j << " = " << pathWeight << " is longer than " << adjMatrix2[i][j] << std::endl;

        }
      }
      // }
    }
  }

//  bool tmp = compareUintMatrices( adjMatrix1, adjMatrix2, numberOfNodes);


  //TODO porownaj matrix z i bez ifa

  std::cout << "A 3 " << std::endl;

  uint32_t steinerWeight = calculateSteiner(
    terminals,
    adjMatrix2,
    allSubsets,
    q);
//  bool tmp = compareUintMatrices( adjMatrix1, adjMatrix2, numberOfNodes);
  for (auto const &x : subsetMap)
  {
    std::cout << "Subset = ";
    std::set<uint32_t> element = x.first;
    for (auto f : element)
    {
      std::cout << f << " ";
    }

    std::cout << std::endl;
    std::cout << "Index = " << x.second << std::endl;
  }
 bool tmp = compareUintMatrices( adjMatrix1, adjMatrix2, numberOfNodes);


 return steinerWeight;
}

//TODO sp oddcycle 3 dobry na naprawe takahasi

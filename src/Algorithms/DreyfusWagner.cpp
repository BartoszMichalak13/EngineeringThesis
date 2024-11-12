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


//TODO clear this file
//TODO check if floyd is maybe faster?
//TODO how about aadjacency matrix instead of list,
//maybe I'll get extra time for no significant increase in memory usage


// #define max 10

// int p[max][max];
// /*
//   All Pairs Shortest Path using Floyd's Algorithm
//  */
// void allpairshort(int a[max][max], int n)
// {
//   int k, i, j;
//   for (k = 0; k < n; k++)
//   {
//     for (i = 0; i < n; i++)
//     {
//       for (j = 0; j < n; j++)
//       {
//         if (a[i][k] + a[k][j] < a[i][j])
//         {
//           a[i][j] = a[i][k] + a[k][j];
//           p[i][j] = k;
//         }
//       }
//     }
//   }
// }

// /*
//  * Storing the shortest path
//  */
// void shortest(int i, int j)
// {
//   int k = p[i][j];
//   if (k > 0)
//   {
//     shortest(i, k);
//     std::cout << "  " << k << "  ";
//     shortest(k, j);
//   }
// }

/*hashTable
Key = {1,2}  Value = 0
Key = {2,3}  Value = 1
{1,4}
{1,2,3}
....
...
{1,2,3,4}
{1,2,3,4} = {1,2}, {2,3}, {3,4}, {1,4}, {2,4}, {1,3} {1,2,3,4}
Target subset size = 3
E = 4
{1,2,3}
{2,3,4}
numberOfSubsets of size r = n!/(n-r)!*r!*/

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
      std::cout << "Adding subset with index" << subsetIndex << "\n";
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

    std::cout << "\n";
    std::cout << "Index = " << x.second << std::endl;
  }

  std::cout<<"Completed subsets";

  return allSubsets;
}










/*
repair
*/
void initializeSteinerTable(std::vector<std::vector<uint32_t>> &steinerDistanceTable, std::vector<std::vector<uint32_t>> &graph, uint32_t numberOfNodes)
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
  {
    for (uint32_t j = 0; j < numberOfNodes; ++j)
    {
      steinerDistanceTable[i][j] = graph[i][j];
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
    // std::vector<std::vector<uint32_t>> graph,
    std::vector<std::set<uint32_t>> allSubsets,
    uint32_t q)
{

  // std::shared_ptr<Graph> self = shared_from_this();

  // Remember that C = Y - {any one element} where Y = {Steiner terminals}
  // Key is index of D and value is distance
  //nodes = {1,2,3,4,5,6,7}
  //D = { {2,3}, {2,4}, {3,4}, {2,3,4} }
  uint32_t sizetable = allSubsets.size() + numberOfNodes; //+ graph.size(); //
  std::vector<std::vector<uint32_t>> steinerDistanceTable(sizetable, std::vector<uint32_t> (numberOfNodes));

  std::cout << "works 1" << std::endl;

  std::vector<std::vector<uint32_t>> adjMatrix = toAdjacencyMatrix();

  std::cout << "works 2" << std::endl;

  //initialise steiner table
  // initializeSteinerTable(steinerDistanceTable, graph, sizetable);
  initializeSteinerTable(steinerDistanceTable, adjMatrix, numberOfNodes);

  uint32_t indexOfSubset = 0;

  std::set<std::set<uint32_t>>::iterator setIterator; // iterator for the "outer" structure

  for (uint32_t m = 2; m < C.size(); ++m)
  {









    // for each subset
    for (uint32_t i = 0; i < allSubsets.size(); ++i)
    {
      std::set<uint32_t> elementOfD = allSubsets[i];
      if(elementOfD.size() == m)
      {
        indexOfSubset = subsetMap[elementOfD];
        // initialise distance from elementOfD to every node in graph
        for (uint32_t k = 0; k < numberOfNodes; ++k)
        {
          steinerDistanceTable[indexOfSubset][k] = std::numeric_limits<uint32_t>::infinity();
        }










        // for every element in graph
        for (uint32_t j = 0; j < numberOfNodes; ++j)
        {
          uint32_t u = std::numeric_limits<uint32_t>::infinity(); ;

          std::set<uint32_t>::iterator nodeIterator;
          for (nodeIterator = elementOfD.begin(); nodeIterator != elementOfD.end(); ++nodeIterator)
          {
            // E is being considered as a single number from the set named elementOfD
            uint32_t E = *nodeIterator;
            uint32_t distanceFromEToJ = steinerDistanceTable[E][j];

            // get index of the subset D-E
            std::cout<<"Getting D-E\n";
            std::cout << "J: "<< j << std::endl;
            std::cout << "Element of D value " << distanceFromEToJ << std::endl;
            std::cout << "Value for E is " << E << std::endl;

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
            uint32_t dist,indexOfSubsetDMinusE;
            if(DMinusE.size() == 1)
            {
              std::set<uint32_t>::iterator node = DMinusE.begin();
              indexOfSubsetDMinusE = *node;
            }
            else
            {
              indexOfSubsetDMinusE = fetchIndexofMapofSets(DMinusE);
            }
            // get the real value here
            std::cout << " indexOfSubsetDMinusE =  " << indexOfSubsetDMinusE << std::endl;
            uint32_t distanceFromDMinusEToJ = steinerDistanceTable[indexOfSubsetDMinusE][j];
            std::cout<<"distanceFromDMinusEToJ = "<< distanceFromDMinusEToJ <<"\n";
            dist = distanceFromEToJ + distanceFromDMinusEToJ;
            std::cout<<"Total Distance = "<< dist <<"\n";
            u = std::min(u, dist);

            std::cout << "U is" << u << std::endl;
          }

          uint32_t val;

          for (uint32_t k = 0; k < numberOfNodes; ++k)
          {
            val = adjMatrix[k][j];
            if (val != 0 )
            {
              steinerDistanceTable[indexOfSubset][k] = std::min(steinerDistanceTable[indexOfSubset][k],val + u);
              std::cout << "I is" << k << std::endl;
              std::cout << "Steiner table data is is" << steinerDistanceTable[indexOfSubset][k] << std::endl;
            }
          }
        }
      }
    }
  }
  uint32_t vm = std::numeric_limits<uint32_t>::infinity();
  uint32_t v = std::numeric_limits<uint32_t>::infinity();

  for (uint32_t j = 0; j < numberOfNodes; j++)
  {
    uint32_t u = std::numeric_limits<uint32_t>::infinity();

    std::set <uint32_t> C1 (C.begin(),C.end());

    for (std::set<uint32_t>::iterator setiter = C1.begin(); setiter != C1.end(); ++setiter)
    {
      uint32_t E = *setiter;
      std::set<uint32_t> subsetCMinusi (C.begin(), C.end());
      uint32_t distEJ = adjMatrix[E][j];
      for (std::set<uint32_t>::iterator iter = subsetCMinusi.begin(); iter != subsetCMinusi.end();)
      {
        if (*iter == E)
        {
          iter = subsetCMinusi.erase(iter);
        }
        else
        {
          ++iter;
        }
      }
      uint32_t indexOfSubsetCMinusi;
      if(subsetCMinusi.size() != 1)
      {
        indexOfSubsetCMinusi = fetchIndexofMapofSets(subsetCMinusi);
      }
      else
      {
        std::set<uint32_t>::iterator node = subsetCMinusi.begin();
        indexOfSubsetCMinusi = *node;
      }

      uint32_t val = steinerDistanceTable[indexOfSubsetCMinusi][j];
      std::cout << " Val: " << val << std::endl;
      std::cout << "E is:  " << E;
      std::cout << "J: "<< j << std::endl;
      std::cout << " distEJ: " << distEJ << std::endl;
      std::cout << " IndexSubset: " << indexOfSubsetCMinusi;

      u = std::min(u,distEJ + val );
      
      std::cout << " U is: " <<u << std::endl;
    }

  std::cout << " U final is:" <<u;
  uint32_t value = adjMatrix[q][j];
  std::cout << "value is: " << value << std::endl;
  if (value < std::numeric_limits<uint32_t>::infinity() && value > 0 )
  {
      v = std::min(v,value+u);
  }
  std::cout << " V is: " << v << std::endl;

  }
  for(uint32_t i = 0 ; i < steinerDistanceTable.size() ; i++ )
  {
    for(uint32_t j = 0; j < numberOfNodes ; j++)
    {
      std::cout << " " << steinerDistanceTable[i][j];
    }
    std::cout << " " << std::endl;
  }


  std::cout << "Minimum Steiner Distance" << v << "\n";
  return v;
}

// void testSubsets()
// {
//   uint32_t n;
//   std::vector<std::set<uint32_t>> allSubsets;
//   std::set<uint32_t>::iterator it;

//   std::cout << "Enter size of the set\n";
//   std::cin >> n;

//   std::vector<uint32_t> arr(n);

//   std::cout << "Enter Elements of the set\n";
//   for (uint32_t i = 0; i < n; i++)
//     std::cin >> arr[i];

//   // allSubsets = generateSubsets(arr, n);

//   for (uint32_t i = 0; i < allSubsets.size(); i++)
//   {
//     std::set<uint32_t> numbersSet = allSubsets[i];
//     for (it = numbersSet.begin(); it != numbersSet.end(); it++)
//     {
//       std::cout << ' ' << *it;
//     }
//     std::cout << "\n";
//   }
// }





















uint32_t Graph::DreyfusWagner(std::vector<uint32_t> terminals) {
  // get graph instance
  std::shared_ptr<Graph> self = shared_from_this();

  std::cout << "HELLO 1" << std::endl;

  printNodeVector(terminals);
  // convert vertices to uint vec
  std::shared_ptr<std::vector<uint32_t>> uintVertices = verticesToUint(vertices, numberOfNodes);
  std::cout << "HELLO 2" << std::endl;

  // calculate all shortest paths
  std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> tmpShortestPaths = AllPairsShortestPath((*uintVertices));
  std::cout << "HELLO 3" << std::endl;


  // int ySize;

  // int a[][10] =  {{0, 3, 1, infi, infi},
  //                 {3, 0 , 7, 5, 1},
  //                 {1, 7 , 0, 2, infi},
  //                 {infi, 5 , 2, 0, 7},
  //                 {infi,1,infi,7,0}
  //                 };
  // allpairshort(a, 5);

  // for (int i = 0; i < 5; i++){        //creating row
  //   graph.push_back(std::vector<int>());
  // }

  // for (int n = 0; n < 5; n++){        //creating columns for the rows
  //   for (int m = 0; m < 5; m++){
  //     graph[m].push_back(0);
  //   }
  // }

////2, 21, 6, 18, 22, 13

  // for (int m = 0; m < 5; m++){        //storing and printing data
  //   for (int n = 0; n < 5; n++){
  //   //  vec[n].push_back(arr[m][n]);
  //     vec[m][n] = a[m][n];
  //   }
  //   std::cout << "\n";
  // }

  // std::cout << "Enter size of terminals\n";
  // std::cin >> ySize;
  // std::vector<int> y(ySize, 0);
  // for (int i = 0; i < ySize; i++)
  // {
  //   std::cin >> y[i];
  // }

  // int q = y[0];
  // y.erase(y.begin());

  uint32_t q = terminals[0];

  std::vector<std::set<uint32_t>> allSubsets = generateSubsets(terminals);

  // uint32_t steinerWeight = calculateSteiner(
  //   terminals,
  //   // graph,
  //   allSubsets,
  //   q);


  /*for(int i = 0 ; i < 5 ; i++ )
  {
      for(int j = 0; j < 5 ; j++)
      {
          cout << " " <<a[i][j];
      }
      cout << " " <<endl;
  }*/

 return 0;
//  return steinerWeight;
}

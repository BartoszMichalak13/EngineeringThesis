#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>
#include <stdint.h>

//TODO segfault na na koniec; do vec kopiuje sie podwojnie, zle kopiuje wierzcholki z wektora, nie wszystkie, a potem problem ze nie znajduje indexu w tablicy / za duzy

Random randomGen;
const uint32_t maxEdgeWeight = 1024;

/*
Used in priorityQueues on Edges
*/
struct EdgeWeightComparatorOnPointers {
  bool operator()(Edge* const& e1, Edge* const& e2) {
    return e1->weight > e2->weight;
  }
};


void printEdgePred(Edge e) {
  if (e.pred == nullptr)
    std::cout << "e.pred doesn't exist; " << std::endl;
  else
    std::cout << "e.pred->end->id " << e.pred->end->id << "; e.pred->start->id " << e.pred->start->id << "; e.pred->weight " << e.pred->weight << "; " << std::endl;
}

void printEdge(Edge e) {
  if (e.end == nullptr) {
    std::cout << "e.end doesn't exist" << std::endl;
    // std::cout << "e.end->id" << e.end->id << std::endl;
    return;
  }
  if (e.start == nullptr) {
    std::cout << "e.start doesn't exist" << std::endl;
    return;
  }

  std::cout << "e.end->id " << e.end->id << "; e.start->id " << e.start->id << "; e.weight " << e.weight << "; ";
  if (e.pred == nullptr)
    std::cout << "e.pred doesn't exist; ";
  else
    std::cout << "e.pred->end->id " << e.pred->end->id << "; e.pred->start->id " << e.pred->start->id << "; e.pred->weight " << e.pred->weight << "; ";
  if (e.succ == nullptr)
    std::cout << "e.succ doesn't exist" << std::endl;
  else
    std::cout << "e.succ->end->id " << e.succ->end->id << "; e.succ->start->id " << e.succ->start->id << "; e.succ->weight" << e.succ->weight << std::endl;
}


void printTmpTreeEdges(Edge ** tmptreeEdges, uint32_t treeSize) {
  for (uint32_t i = 0; i < treeSize; ++i) {
    printEdge((*tmptreeEdges)[i]);
  }
}

/*
returns index from vector which holds value
*/
int32_t findInVector(uint32_t value, std::vector<uint32_t> vec) {
  for (uint32_t i = 0; i < vec.size(); ++i)
    if (vec.at(i) == value)
      return i;
  return -1;
}

/*
returns idex from array of Nodes which holds value
*/
int32_t findInArray(uint32_t value, Node * array, uint32_t arraySize) {
  for (uint32_t i = 0; i < arraySize; ++i)
    if (array[i].id == value)
      return i;
  return -1;
}

/*
returns idex from array of Nodes which holds value
*/
int32_t findInTmpTree(uint32_t edgeStart, uint32_t edgeEnd, Edge* * array, uint32_t arraySize) {
  for (uint32_t i = 0; i < arraySize; ++i)
    if (array[i]->start->id == edgeStart && array[i]->end->id == edgeEnd)
      return i;
  return -1;
}

/*
finds an edge in adjacencyList and returns pointer to it
*/
Edge* findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<Edge*>* adjacencyList) {
  for(uint32_t i = 0; i < adjacencyList[edgeStart].size(); ++i)
    if (adjacencyList[edgeStart].at(i)->end->id == edgeEnd)
      return adjacencyList[edgeStart].at(i);
  return nullptr;
}

// Edge* findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<Edge> edgeVec) {
//   for(uint32_t i = 0; i < edgeVec.size(); ++i)
//     if ((edgeVec.at(i).end->id == edgeEnd && edgeVec.at(i).start->id == edgeStart) ||
//         (edgeVec.at(i).start->id == edgeEnd && edgeVec.at(i).end->id == edgeStart))
//       return &edgeVec.at(i);
//   return nullptr;
// }

Graph::Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<Edge*>[numberOfNodes];
  this->vertices = new Node[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    vertices[i] = Node(i);
}


Graph::Graph(Edge * edges, Node* * nodes, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  std::cout << "numberOfNodes: " << numberOfNodes << std::endl;
  std::cout << "numberOfEdges: " << numberOfEdges << std::endl;
  // std::cout << "nodes: " << sizeof(nodes) << std::endl;
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<Edge*>[numberOfNodes];
  this->vertices = new Node[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    this->vertices[i] = Node(nodes[i]->id);
  for (uint32_t i = 0; i < numberOfEdges; ++i) {
    std::cout << "a " << i << std::endl;
    int32_t idx1 = findInArray(edges[i].start->id, vertices, numberOfNodes);
    if (idx1 < 0)
      std::cout << "Error node not found in Constructor: " << edges[i].start->id << std::endl;
    std::cout << "idx1 " << idx1 << std::endl;
    int32_t idx2 = findInArray(edges[i].end->id, vertices, numberOfNodes);
    if (idx2 < 0)
      std::cout << "Error node not found in Constructor: " << edges[i].end->id << std::endl;
    std::cout << "idx2 " << idx2 << std::endl;
    Node * start = &vertices[idx1];
    std::cout << "d " << i << std::endl;
    Node * end = &vertices[idx2];
    std::cout << "e " << i << std::endl;
    Edge e1(start, edges[i].weight, end);
    std::cout << "f " << i << std::endl;
    Edge e2(end, edges[i].weight, start);
    std::cout << "g " << i << std::endl;
    this->adjacencyList[idx1].push_back(&e1);
    std::cout << "h " << i << std::endl;
    this->adjacencyList[idx2].push_back(&e2);
    std::cout << "i " << i << std::endl;
    printEdge(e1);
    printEdge(e2);
    std::cout << "numberOfPushedEdges: " << i*2 << std::endl;
  }
}

// Graph::Graph(Edge * edges, Node * nodes, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
//   std::cout << "numberOfNodes: " << numberOfNodes << std::endl;
//   std::cout << "numberOfEdges: " << numberOfEdges << std::endl;
//   // std::cout << "nodes: " << sizeof(nodes) << std::endl;
//   this->printFlag = printFlag;
//   this->numberOfNodes = numberOfNodes;
//   this->numberOfEdges = numberOfEdges;
//   this->adjacencyList = new std::vector<Edge>[numberOfNodes];
//   this->vertices = new Node[numberOfNodes];
//   for (uint32_t i = 0; i < numberOfNodes; ++i)
//     this->vertices[i] = Node(nodes[i].id);
//   for (uint32_t i = 0; i < numberOfEdges; ++i) {
//     std::cout << "a " << i << std::endl;
//     int32_t idx1 = findInArray(edges[i].start->id, vertices, numberOfNodes);
//     if (idx1 < 0)
//       std::cout << "Error node not found in Constructor: " << edges[i].start->id << std::endl;
//     std::cout << "idx1 " << idx1 << std::endl;
//     int32_t idx2 = findInArray(edges[i].end->id, vertices, numberOfNodes);
//     if (idx2 < 0)
//       std::cout << "Error node not found in Constructor: " << edges[i].end->id << std::endl;
//     std::cout << "idx2 " << idx2 << std::endl;
//     Node * start = &vertices[idx1];
//     std::cout << "d " << i << std::endl;
//     Node * end = &vertices[idx2];
//     std::cout << "e " << i << std::endl;
//     Edge e1(start, edges[i].weight, end);
//     std::cout << "f " << i << std::endl;
//     Edge e2(end, edges[i].weight, start);
//     std::cout << "g " << i << std::endl;
//     this->adjacencyList[idx1].push_back(e1);
//     std::cout << "h " << i << std::endl;
//     this->adjacencyList[idx2].push_back(e2);
//     std::cout << "i " << i << std::endl;
//     printEdge(e1);
//     printEdge(e2);
//     std::cout << "numberOfPushedEdges: " << i*2 << std::endl;
//   }
// }

Graph::~Graph() {
  delete[] adjacencyList;
  // adjacencyList = nullptr;
  delete[] vertices;
  // vertices = nullptr;
}

Graph dummyGraph() {
  return Graph(0,0,0);
}

void Graph::addEdge(uint32_t node1Id, uint32_t node2Id) {
  // We choose them, so that they're smaller than numberOfNodes
  const uint32_t weight = randomGen.generateRandomNumber(1, maxEdgeWeight);
  if(weight == 0)
    std::cout << "ERRRRRRRRRRROR\n";
  Edge e1(&vertices[node1Id], weight, &vertices[node2Id]);
  Edge e2(&vertices[node2Id], weight, &vertices[node1Id]);
  this->adjacencyList[node1Id].push_back(&e1);
  this->adjacencyList[node2Id].push_back(&e2);
}

void printAdajcencyList(std::vector<Edge*>* adjList, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    std::cout << "node "<< i << std::endl;
    std::cout << "\t";
    // for (const Edge& e: adjList[i]) {
    for (uint32_t j = 0; j < adjList[i].size(); ++j) {
      Edge * e = adjList[i].at(j);
      std::cout << e->end->id << " cost(" << e->weight << ")";
      if (e->weight == 0)
        std::cout << "*** ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
}

/*
Lists should have the same number of vectors.
@return returns true if they are different
*/
bool compareAdajcencyLists(std::vector<Edge*>* adjList1, std::vector<Edge*>* adjList2, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    if (adjList1[i].size() != adjList2[i].size()) {
      std::cout << "Different number of elementsa at "<< i << std::endl;
      for (uint32_t j = 0; j < adjList1[i].size(); ++j)
        std::cout << adjList1[i].at(j)->end->id << ",";
      std::cout << std::endl;
      for (uint32_t j = 0; j < adjList2[i].size(); ++j)
        std::cout << adjList2[i].at(j)->end->id << ",";
      std::cout << std::endl;
      return true;
    } else {
      for (uint32_t j = 0; j < adjList1[i].size(); ++j) {
        if (adjList1[i].at(j)->end != adjList2[i].at(j)->end && adjList1[i].at(j)->start != adjList2[i].at(j)->start)
        {
          std::cout << "Difference in row " << i << " at element "<< j << " of value: "
            << adjList1[i].at(j)->start->id << ","<< adjList1[i].at(j)->end->id << " vs "
            << adjList2[i].at(j)->start->id << ","<< adjList2[i].at(j)->end->id  << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}

/*
Generates vector of idx's of terminals.
*/
std::vector<uint32_t> Graph::generateTerminals(uint32_t numberOfTerminals) {
  std::vector<uint32_t> terminals;
  uint32_t term = randomGen.generateRandomNumber(0, numberOfNodes - 1);
  terminals.push_back(term);
  for (uint32_t i = 1; i < numberOfTerminals; ++i) {
    while (findInVector(term, terminals) > -1) {
      term = randomGen.generateRandomNumber(0, numberOfNodes - 1);
    }
    terminals.push_back(term);
  }
  return terminals;
}

void Graph::bfs() {
  std::queue<Node*> toVisit;
  Node * currentNode = &vertices[0];
  toVisit.push(currentNode);
  while (!toVisit.empty()) {
    currentNode = toVisit.front();
    toVisit.pop();
    currentNode->visited = true;
    for (const auto& edge : adjacencyList[currentNode->id]) {
      if(!edge->end->visited) {
        toVisit.push(edge->end);
        edge->end->visited = true;
      }
      //printVisitedStatus();
    }
  }
} 

/*
copies adjlist of a graph to given copyOfAdjacencyList
*/
void Graph::copyAdjacencyList(std::vector<Edge*>* copyOfAdjacencyList) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    // for (const Edge& edge : adjacencyList[i]) {
    for (uint32_t j = 0; j < adjacencyList[i].size(); ++j) {
      Edge * edge = adjacencyList[i].at(j);
      Node* start = &vertices[edge->start->id];
      Node* end = &vertices[edge->end->id];
      Edge copyOfEdge(start, edge->weight, end, edge->pred, edge->succ);
      copyOfAdjacencyList[i].push_back(&copyOfEdge);
    }
}

//TODO: make sure generator, generates in range 1-something aka non zero
/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void resetCopyOfAdjacencyList(std::vector<Edge*>*& localCopyOfAdjacencyList, Graph * g) {
  for (uint32_t i = 0; i < g->numberOfNodes; ++i) {
    for (uint32_t j = 0; j < localCopyOfAdjacencyList[i].size(); ++j) {
      if (localCopyOfAdjacencyList[i].at(j)->weight)
      // { // if non zero edge - reset it to default (adjacencyList value)
        // Edge * pred = localCopyOfAdjacencyList[i].at(j).pred;
        // Edge * succ = localCopyOfAdjacencyList[i].at(j).succ;
        localCopyOfAdjacencyList[i].at(j)->weight = g->adjacencyList[i].at(j)->weight; //TODO: does it work as intended?
        // localCopyOfAdjacencyList[i].at(j).pred = pred;
        // localCopyOfAdjacencyList[i].at(j).succ = succ;
      // } else { // weight == 0 aka edge in tree
      //   // Edge * pred = &localCopyOfAdjacencyList[i].at(j);
      //   // localCopyOfAdjacencyList[i].at(j).pred = pred;//nullptr // TODO: is this correct way to do this?
      //   // localCopyOfAdjacencyList[i].at(j).succ = nullptr;
      // }
      localCopyOfAdjacencyList[i].at(j)->end->visited = false;
      localCopyOfAdjacencyList[i].at(j)->start->visited = false;
    }
  }
}

/*
TODO: CHECK IF WORKS CORRECTLY
Return array of Edges with non empty fields passed to second parameter.
*/
void cropArraytoVector(Edge* * tmptreeEdges, std::vector<Edge>& vec, uint32_t treeSize) {
  for (uint32_t i = 0; i < treeSize; ++i)
    vec.push_back((*tmptreeEdges)[i]);
}

//TODO two functions below: check if it should be as pointer or as a value

/*
TODO: CHECK IF WORKS CORRECTLY
Return array of Edges with non empty fields passed to second parameter.
*/
void cropArray(Edge* * tmptreeEdges, Edge* *treeEdges, uint32_t treeSize) {
  for (uint32_t i = 0; i < treeSize; ++i)
    *(treeEdges)[i] = *(tmptreeEdges)[i];
}
/*
Returns array of Nodes from array of Edges (pass number of edges as treeSize, function knows that |V| = |E| + 1)
*/
// Node * getNodesFromEdges(std::vector<Edge> treeEdges, uint32_t treeSize) {
//   Node * nodeArray = new Node[treeSize + 1];
//   nodeArray[0] = *treeEdges.at(0).start;
//   for (uint32_t i = 0; i < treeSize; ++i)
//     nodeArray[i + 1] = *treeEdges.at(i).end;
//   return nodeArray;
// }

/*
Returns array of Nodes from array of Edges. treeSize is number of edges
*/
uint32_t getNodesFromEdges(std::vector<Edge> treeEdges, Node ** &treeNodes, uint32_t treeSize) {
  uint32_t currentNumberOfNodes = 0;
  for (uint32_t i = 0; i < treeSize; ++i) {
    if (findInArray(treeEdges.at(i).start->id, *treeNodes, currentNumberOfNodes) == -1) {
      treeNodes[currentNumberOfNodes] = treeEdges.at(i).start;
      std::cout << treeNodes[currentNumberOfNodes]->id << ", ";
      ++currentNumberOfNodes;
    }
    if (findInArray(treeEdges.at(i).end->id, *treeNodes, currentNumberOfNodes) == -1) {
      treeNodes[currentNumberOfNodes] = treeEdges.at(i).end;
      std::cout << treeNodes[currentNumberOfNodes]->id << ", ";
      ++currentNumberOfNodes;
    }
  }
  std::cout << std::endl;
  return currentNumberOfNodes;
}

/*
Finds both instances of edge and updates their weights.
*/
void updateWeight(Edge * e, uint32_t weight, std::vector<Edge*>* adjList) {
  e->weight = weight;
  Edge * edge = findEdge(e->end->id, e->start->id, adjList);
  if (edge == nullptr) {
    std::cout << "Couldn't find edge: "; printEdge(*e);
  } else {
    edge->weight = weight;
  }
}

/*
Finds both instances of edge and updates their pred.
*/
void updatePred(Edge * e, Edge * pred, std::vector<Edge*>* adjList) {
  e->pred = pred;
  findEdge(e->end->id, e->start->id, adjList)->pred = pred;
}

/*
Finds both instances of edge and creates pred loop.
*/
void updatePredToLoop(Edge * e, std::vector<Edge*>* adjList) {
  e->pred = e;
  Edge * e2 = findEdge(e->end->id, e->start->id, adjList);
  e2->pred = e2;
}

/*
Finds both instances of edge and updates their succ.
*/
void updateSucc(Edge * e, Edge * succ, std::vector<Edge*>* adjList) {
  e->succ = succ;
  findEdge(e->end->id, e->start->id, adjList)->succ = succ;
}

/*difference between prim - 
1. we update weights in edges, based on distance to beg. node 
2. in this variation of Dijkstra we search just untill we meet node we search for
3. Terminals should be narrowed down after each found one
KIEDY ZNAJDE WAGI W NOWYM GRAFIE NA 0000000000000000000

usuwam wiekszosc pointerow, gdyz tworze nowe drzewo, a nie chce uszkadzac starego

modify edges, give pointer to pred

@param terminals should be of positive size
*/
Graph Graph::Dijkstra(std::vector<uint32_t> terminals) {// do it with priority Queue, maybe my own immplementation?
// Original terminals: 27, 0, 35, 20, 14, 6, 28, 34, 8, 25,
  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  std::cout << "Original terminals: ";
  for (uint32_t i = 0; i < originalTerminals.size(); ++i)
    std::cout << originalTerminals.at(i) << ", ";
  std::cout << std::endl;

  //init
  resetVisitedStatus();
  std::vector<Edge*>* localCopyOfAdjacencyList = new std::vector<Edge*>[numberOfNodes];
  copyAdjacencyList(localCopyOfAdjacencyList);

  printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);

  Edge * tmptreeEdges[numberOfNodes - 1]; // no chyba nie tyle I guess ale maxymalnie tyle

  uint32_t treeSize = 0;
  bool foundTerminal = false;
  std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparatorOnPointers> toVisit;

  // for (const Edge& edge : localCopyOfAdjacencyList[terminals.at(0)]) {
  for (uint32_t j = 0; j < localCopyOfAdjacencyList[0].size(); ++j) {
    Edge * edge = localCopyOfAdjacencyList[0].at(j);
    std::cout << "terminals.at(0) " << terminals.at(0) << std::endl;
    printEdge(edge);
    toVisit.push(edge);
  }

  vertices[terminals.at(0)].visited = true;
  terminals.erase(terminals.begin());
  do { //TMP NOTE - should be 9
    std::cout << "dbg9 " << std::endl;
    if (compareAdajcencyLists(adjacencyList, localCopyOfAdjacencyList, numberOfNodes)) {
      std::cout << "Comparing Adajcency Lists Failed" << std::endl;
      return dummyGraph();
    }
    std::cout << "dbg10 " << std::endl;

    if (foundTerminal) {
      while(!toVisit.empty())
        toVisit.pop();
      std::cout << "dbg11 " << std::endl;
      printTmpTreeEdges(tmptreeEdges, treeSize);
      std::cout << "dbg11.25 " << std::endl;
      std::vector<Edge> tmpEdges;
      cropArraytoVector(tmptreeEdges, tmpEdges, treeSize);
      Node* * tmpNodes = new Node*[treeSize + 1];
      for (uint32_t i = 0; i < tmpEdges.size(); ++i) {
        printEdge(tmpEdges.at(i));
      }
      std::cout << "dbg11.5 " << treeSize << std::endl;

      uint32_t currentNumberOfNodes = getNodesFromEdges(tmpEdges, tmpNodes, treeSize);
      std::cout << "currentNumberOfNodes " << currentNumberOfNodes << std::endl;
      Edge* * e = new Edge*[treeSize + 1];
      std::cout << "debug tsize " << treeSize << std::endl;
      printAdajcencyList(localCopyOfAdjacencyList, numberOfNodes);
      for (uint32_t i = 0; i < treeSize + 1; ++i) {//for each node in current tree
        // for (const Edge& edge : localCopyOfAdjacencyList[tmpNodes[i]->id]) {
        for (uint32_t j = 0; j < localCopyOfAdjacencyList[tmpNodes[i]->id].size(); ++j) {
          Edge * edge = localCopyOfAdjacencyList[tmpNodes[i]->id].at(j);
          if (!edge->weight){
            e[i] = edge;
            std::cout << "e[" << i << "]->weight " << e[i]->weight << std::endl;
          }
        }
        std::cout << "debug " <<  i << std::endl;
        std::cout << "tmpNodes[" << i << "]->id " << tmpNodes[i]->id << std::endl;
        // add all other edges and set thier pred as edge from tree
        // for (Edge& edge : localCopyOfAdjacencyList[tmpNodes[i]->id]) {
        for (uint32_t j = 0; j < localCopyOfAdjacencyList[tmpNodes[i]->id].size(); ++j) {
          Edge * edge = localCopyOfAdjacencyList[tmpNodes[i]->id].at(j);
          if (edge->weight) // normal procedure
            updatePred(edge, e[i], localCopyOfAdjacencyList);
          else // if it's 0 it's part of a tree
            updatePredToLoop(edge, localCopyOfAdjacencyList);

          updateWeight(edge, (edge->weight + e[i]->weight), localCopyOfAdjacencyList);
          toVisit.push(edge);
        }
      }
      std::cout << "dbg12 " << std::endl;

      delete[] tmpNodes;
      // tmpNodes = nullptr;
      delete[] e;
      // e = nullptr;

      foundTerminal = false;
    }

    //main loop
    while(!foundTerminal) { // O(n^2)
      //remoeve all edges that lead to already visited nodes
      if (!toVisit.empty())
        while (toVisit.top()->end->visited)
          toVisit.pop();

      if (toVisit.empty()) {
        std::cout << "End in loop Dijkstra; dummy graph" << std::endl;
        return dummyGraph();
      }

      Edge * e = toVisit.top();
      toVisit.pop();
      e->end->visited = true;
      uint32_t nextNodeIndex = e->end->id;

      // check if we have found terminal
      int32_t idx = findInVector(nextNodeIndex, terminals);
      std::cout << "WTF idx:" << idx << std::endl;
      printEdge(*e);
      if (idx > -1) {
        foundTerminal = true;
        //remove found terminal from the list
        terminals.erase(terminals.begin() + idx);
        std::cout << "WTF " << std::endl;
        printEdge(*e);

        //TODO - wouldnt this work if I just add the edge, and it;s pred in loop?

        /*
        AT THIS POINT
        1. zero the path
        2.  a) reset all weights != 0
            b) reset all visited status
            c) reset all pred
            d) reset all succc
        3. Reset Queue
        4. Rerun algorithm with new local adjacencyList (changed in 1. and 2.)
        */
        if (e->pred == nullptr) {
          std::cout << "dbg1 " << std::endl;
          // std::cout << "e->end->id " << e->end->id << std::endl;
          // std::cout << "e->start->id " << e->start->id << std::endl;
          printEdge(*findEdge(e->start->id, e->end->id, adjacencyList));
          Edge * edg = findEdge(e->start->id, e->end->id, adjacencyList);
          if(findInTmpTree(edg->start->id, edg->end->id, tmptreeEdges, treeSize) == -1) {
            std::cout << "added " << std::endl;
            printEdge(*edg);
            tmptreeEdges[treeSize] = edg;
            ++treeSize;
            printTmpTreeEdges(tmptreeEdges, treeSize);
          }
          std::cout << "dbg2 " << std::endl;
          if (tmptreeEdges[treeSize] == nullptr) {
            std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
            return dummyGraph();
          }
          updateWeight(e, 0, localCopyOfAdjacencyList);
          updatePredToLoop(e, localCopyOfAdjacencyList);
        } else {
          // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
          while (e->pred->start->id != e->start->id) {
            std::cout << "dbg1 " << std::endl;
            std::cout << "HELLO " << std::endl;
            std::cout << (e->pred != nullptr) << std::endl;
            std::cout << (e->pred->start->id != e->start->id) << std::endl;
            // std::cout << "e->end->id " << e->end->id << std::endl;
            // std::cout << "e->start->id " << e->start->id << std::endl;
            std::cout << "dbg1.5 " << std::endl;
            updateWeight(e, 0, localCopyOfAdjacencyList);
            printEdge(*e);
            Edge * oldPred = e->pred;
            updatePredToLoop(e, localCopyOfAdjacencyList);
            printEdge(*findEdge(e->start->id, e->end->id, adjacencyList));
            Edge * edg = findEdge(e->start->id, e->end->id, adjacencyList);
            if(findInTmpTree(edg->start->id, edg->end->id, tmptreeEdges, treeSize) == -1) {
              std::cout << "added " << std::endl;
              printEdge(*edg);
              tmptreeEdges[treeSize] = edg;
              ++treeSize;
              printTmpTreeEdges(tmptreeEdges, treeSize);
            }
            std::cout << "dbg2 " << std::endl;
              printEdge(*e);

            if (tmptreeEdges[treeSize] == nullptr) {
              std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
              return dummyGraph();
            }
            std::cout << "dbg3 " << std::endl;
            std::cout << "dbg4 " << std::endl;
            std::cout << "dbg5 " << std::endl;
            printEdge(*e);
            if (oldPred == nullptr) {
              std::cerr << "oldPred == nullptr" << std::endl;
              break;
            } else {
              std::cout << "dbg5.1 " << std::endl;
              printEdge(*e);
              std::cout << "dbg5.2 " << std::endl;
              printEdgePred(*e);
              std::cout << "dbg5.2 " << std::endl;
              e = oldPred;
              if (e->pred == nullptr){
                std::cout << "e->pred == nullptr in while "<< std::endl;
                break;
              }
              printEdge(*e);
              printEdgePred(*e);
            }
            std::cout << "dbg6 " << std::endl;
          } // untill we reach beginning of the path
          printEdge(*e);
          updateWeight(e, 0, localCopyOfAdjacencyList);
          updatePredToLoop(e, localCopyOfAdjacencyList);
          Edge * edg = findEdge(e->start->id, e->end->id, adjacencyList);
          if(findInTmpTree(edg->start->id, edg->end->id, tmptreeEdges, treeSize) == -1) {
            std::cout << "added " << std::endl;
            printEdge(*edg);
            tmptreeEdges[treeSize] = edg;
            ++treeSize;
            printTmpTreeEdges(tmptreeEdges, treeSize);
          }
          if (tmptreeEdges[treeSize] == nullptr) {
            std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
            return dummyGraph();
          }
        }
        std::cout << "dbg6.5 " << std::endl;
// Original terminals: x38, 19, 9, 28, 13, x7, x11, 14, 34, 33,

        std::cout << "dbg7 " << std::endl;
        resetCopyOfAdjacencyList(localCopyOfAdjacencyList, this);
        std::cout << "dbg8 " << std::endl;
        break;
      }

      for (uint32_t j = 0; j < localCopyOfAdjacencyList[nextNodeIndex].size(); ++j) {
        Edge * edge = localCopyOfAdjacencyList[nextNodeIndex].at(j);
      // for (Edge& edge : localCopyOfAdjacencyList[nextNodeIndex]) {
        if (!edge->end->visited) {
          if (edge->weight) // normal procedure
            updatePred(edge, e, localCopyOfAdjacencyList);
          else // if it's 0 it's part of a tree
            updatePredToLoop(edge, localCopyOfAdjacencyList);
          updateWeight(edge, (edge->weight + e->weight), localCopyOfAdjacencyList);
          toVisit.push(edge);
        }
      }
    }
  } while(!terminals.empty()); // O(k)

  uint64_t totalWeight = 0;
  for (uint32_t i = 0; i < treeSize; ++i)
    totalWeight += tmptreeEdges[i]->weight;
  //TODO: print edges in separate function
  // print edges
  if (printFlag && numberOfNodes < 1000) {
    for (uint32_t i = 0; i < treeSize; ++i)
      std::cout << tmptreeEdges[i]->start->id << "->" <<  tmptreeEdges[i]->end->id << "; ";
    std::cout << std::endl;
    std::cout << "Original terminals: ";
    for (uint32_t i = 0; i < originalTerminals.size(); ++i)
      std::cout << originalTerminals.at(i) << ", ";
    std::cout << std::endl;
  }
  std::cout << "Dijkstra totalWeight = " << totalWeight << std::endl;


  std::cout << "Edges in  tmptreeEdges: " << std::endl;
  for (uint32_t i = 0; i < treeSize; ++i) {
    printEdge(*tmptreeEdges[i]);
  }

  std::vector<Edge> treeEdgesVec;
  cropArraytoVector(tmptreeEdges, treeEdgesVec, treeSize);

  std::cout << "Edges in  treeEdgesVec: " << std::endl;
  for (uint32_t i = 0; i < treeEdgesVec.size(); ++i) {
    printEdge(treeEdgesVec.at(i));
  }

  Edge * treeEdges = new Edge[treeSize];
  for (uint32_t i = 0; i < treeEdgesVec.size(); ++i) {
    treeEdges[i] = treeEdgesVec.at(i);
  }
  for (uint32_t i = 0; i < treeEdgesVec.size(); ++i) {
    printEdge(treeEdges[i]);
  }
  // cropArray(tmptreeEdges, &treeEdges, treeSize);
  Node* * treeNodes = new Node*[treeSize + 1];
  uint32_t currentNumberOfNodes = getNodesFromEdges(treeEdgesVec, treeNodes, treeSize);
    std::cout << "currentNumberOfNodes "<< currentNumberOfNodes << std::endl;
    std::cout << "dupa 123" << std::endl;

  for (uint32_t i = 0; i < treeSize + 1; ++i) {
    std::cout << treeNodes[i]->id << ", ";
  }
  std::cout << std::endl;

  for (uint32_t i = 0; i < originalTerminals.size(); ++i) {
    bool found = false;
    for (uint32_t j = 0; j < treeSize+1; ++j) {
      if (originalTerminals.at(i) == treeNodes[j]->id) {
          std::cout << "found: " << originalTerminals.at(i) << std::endl;
        found = true;
        break;
      }
    }
    if (!found) {
      std::cout << std::endl;
      std::cout << "missing: " << originalTerminals.at(i) << std::endl;
    }
  }
  std::cout << std::endl;


  // getNodesFromEdges(treeEdges, &treeNodes, treeSize);
  Graph steinerTree = Graph(treeEdges, treeNodes, static_cast<uint32_t>(treeSize + 1), treeSize, printFlag);
  std::cout << "dupa 321" << std::endl;

  delete[] localCopyOfAdjacencyList;
  // localCopyOfAdjacencyList = nullptr;
  std::cout << "dupa 321" << std::endl;

  delete[] treeEdges;
  // treeEdges = nullptr;
  std::cout << "dupa 321" << std::endl;

  delete[] treeNodes;
  // treeNodes = nullptr;
  std::cout << "dupa 321" << std::endl;

  return steinerTree;
}

/*
Initialize a tree with a single vertex, chosen arbitrarily from the graph.
Grow the tree by one edge: Of the edges that connect the tree to vertices not yet in the tree, find the minimum-weight edge, and transfer it to the tree.
Repeat step 2 (until all vertices are in the tree).
*/
Graph Graph::PrimMST() // do it with priority Queue, maybe my own immplementation?
{
  resetVisitedStatus();
  std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparatorOnPointers> toVisit;
  // for (Edge& edge : adjacencyList[0]) //init
  for (uint32_t j = 0; j < adjacencyList[0].size(); ++j) {
    Edge * edge = adjacencyList[0].at(j);
    toVisit.push(edge);
  }
  vertices[0].visited = true;
  
  Edge * tmptreeEdges[numberOfNodes - 1];
  uint64_t totalWeight = 0;
  for (uint32_t treeSize = 0; treeSize < numberOfNodes - 1; ++treeSize) {
    while (toVisit.top()->end->visited) //remoeve all edges that lead to already visited nodes
      toVisit.pop();

    if (toVisit.size() == 0) {
      std::cout << "End in loop PimMST; dummy graph" << std::endl;
      return dummyGraph();
    }
    Edge * e = toVisit.top();
    tmptreeEdges[treeSize] = e;
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    for (uint32_t j = 0; j < adjacencyList[nextNodeIndex].size(); ++j) {
      Edge * edge = adjacencyList[nextNodeIndex].at(j);
      if (!edge->end->visited) {
        toVisit.push(edge);
      }
    }
    totalWeight += e->weight;
  }
  if (printFlag && numberOfNodes < 1000) {
    for (uint32_t i = 0; i < numberOfNodes - 1; ++i)
      std::cout << tmptreeEdges[i]->start->id << "->" <<  tmptreeEdges[i]->end->id << "; ";
    std::cout << std::endl;
  }
  std::cout << "PimMST totalWeight = " << totalWeight << std::endl;

  std::cout << "Dummy graph primMST" << std::endl;
  return dummyGraph();
}



void Graph::printData() {
  for (uint32_t i = 0; i < numberOfNodes; ++i) {
    std::cout << "node " << i << " is connected to:\n\t"; 
    // for (const auto& edge : adjacencyList[i]) {
    for (uint32_t j = 0; j < adjacencyList[i].size(); ++j) {
      Edge * edge = adjacencyList[i].at(j);
      std::cout << static_cast<uint32_t>(edge->end->id) << " cost(" << edge->weight << ") ";
      if(edge->weight == 0)
        std::cout << "ERROR;";
    }
    std::cout << std::endl;
  }
}

bool Graph::isConnected()
{
  bfs();
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    if (!vertices[i].visited)
      return false;
  return true;
}

void Graph::resetVisitedStatus()
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    vertices[i].visited = false;
}

void Graph::printVisitedStatus()
{
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    std::cout << "node " << i << " visited = " << vertices[i].visited << std::endl;
}

bool Graph::checkIfEdgeExists(uint32_t node1Id, uint32_t node2Id) 
{
  // for (const auto& edge : adjacencyList[node1Id])

  for (uint32_t j = 0; j < adjacencyList[node1Id].size(); ++j) {
    Edge * edge = adjacencyList[node1Id].at(j);
    if (edge->end->id == node2Id)
      return true;
  }
  return false;
}

void generateGraph(uint32_t numberOfNodes, uint32_t numberOfEdges, float density, bool  printFlag)
{
  Graph graph = Graph(numberOfNodes, numberOfEdges, printFlag);
  const uint32_t range = numberOfNodes - 1;

  for (uint32_t i = 0; i < numberOfEdges; ++i)
  {
    bool cannotMakeEdge = true;
    while(cannotMakeEdge) 
    {
      uint32_t node1Id = randomGen.generateRandomNumber(0, range);
      uint32_t node2Id = randomGen.generateRandomNumber(0, range);
      cannotMakeEdge = graph.checkIfEdgeExists(node1Id, node2Id);
      if (!cannotMakeEdge && node1Id != node2Id) //no multi edges, no self-loops
        graph.addEdge(node1Id, node2Id);
    }
  }
  if (printFlag && numberOfNodes < 1000)
  {
    graph.printData();
  }
  if(graph.isConnected()) {
    std::cout << "Graph is connected" << std::endl;
  } else {
    std::cout << "Graph is NOT connected" << std::endl;
    //TMP SOLUTION
    return;
  }
  graph.resetVisitedStatus();
  Graph mst = graph.PrimMST();
  graph.resetVisitedStatus();

  //TODO: maybe make them not random?
  //TMP SOLUTION
  uint32_t numberOfTerminals = std::round(numberOfNodes / 4);
  std::cout << "dupa" << std::endl;
  std::vector<uint32_t> terminals = graph.generateTerminals(numberOfTerminals);
  std::cout << "dupa" << std::endl;

  Graph steinerTree = graph.Dijkstra(terminals);
  std::cout << "dupa" << std::endl;
  std::cout << "steinerTree.numberOfNodes " << steinerTree.numberOfNodes << std::endl;
  std::cout << "dupa" << std::endl;
  std::cout << "steinerTree.adjacencyList[0].size() " << steinerTree.adjacencyList[0].size() << std::endl;
  std::cout << "graph.adjacencyList[0].size() " << graph.adjacencyList[0].size() << std::endl;
  std::cout << "dupa" << std::endl;
  for (uint32_t i = 0; i < steinerTree.numberOfNodes; ++i)
    for (uint32_t j = 0; j < steinerTree.adjacencyList[i].size(); ++j) {
      printEdge(steinerTree.adjacencyList[i].at(j));
      printEdgePred(steinerTree.adjacencyList[i].at(j));
    }
  std::cout << "dupa" << std::endl;
}



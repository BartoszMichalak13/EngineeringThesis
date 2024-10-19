#include "graph.hpp"
#include "random.hpp"
#include <iostream>
#include <queue>
#include <functional>
#include <stdint.h>

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
finds an edge in adjacencyList and returns pointer to it
*/
Edge* findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<Edge>* adjacencyList) {
  for(uint32_t i = 0; i < adjacencyList[edgeStart].size(); ++i)
    if (adjacencyList[edgeStart].at(i).end->id == edgeEnd)
      return &adjacencyList[edgeStart].at(i);
  return nullptr;
}

Edge* findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<Edge> edgeVec) {
  for(uint32_t i = 0; i < edgeVec.size(); ++i)
    if ((edgeVec.at(i).end->id == edgeEnd && edgeVec.at(i).start->id == edgeStart) ||
        (edgeVec.at(i).start->id == edgeEnd && edgeVec.at(i).end->id == edgeStart))
      return &edgeVec.at(i);
  return nullptr;
}

Graph::Graph(uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<Edge>[numberOfNodes];
  this->vertices = new Node[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    vertices[i] = Node(i);
}
//TODO change EDGE constructor so its less confusing
Graph::Graph(Edge * edges, Node * nodes, uint32_t numberOfNodes, uint32_t numberOfEdges, bool printFlag) {
  this->printFlag = printFlag;
  this->numberOfNodes = numberOfNodes;
  this->numberOfEdges = numberOfEdges;
  this->adjacencyList = new std::vector<Edge>[numberOfNodes];
  this->vertices = new Node[numberOfNodes];
  for (uint32_t i = 0; i < numberOfNodes; ++i)
    this->vertices[i] = Node(nodes[i].id);

  // std::vector<Edge> edgesVec;
  for (uint32_t i = 0; i < numberOfEdges; ++i) {
    // Node * start = &vertices[findInArray(edges[i].start->id, vertices, numberOfNodes)];
    // Node * end = &vertices[findInArray(edges[i].end->id, vertices, numberOfNodes)];
    // Edge e(end, edges[i].weight, start, edges[i].pred, edges[i].succ);
    // edgesVec.push_back(e);

    Node * start = &vertices[findInArray(edges[i].start->id, vertices, numberOfNodes)];
    Node * end = &vertices[findInArray(edges[i].end->id, vertices, numberOfNodes)];
    // Edge * pred = findEdge(edges[i].pred->start->id,edges[i].pred->end->id,edgesVec);

    Edge e1(end, edges[i].weight, start, edges[i].pred, edges[i].succ);
    Edge e2(start, edges[i].weight, end, edges[i].pred, edges[i].succ);
    uint32_t idx1 = findInArray(edges[i].start->id, vertices, numberOfNodes);
    uint32_t idx2 = findInArray(edges[i].end->id, vertices, numberOfNodes);
    this->adjacencyList[idx1].push_back(e1);
    this->adjacencyList[idx2].push_back(e2);
  }
  //TODO REPAIR PREDECESSORS
  // for (uint32_t i = 0; i < numberOfEdges; ++i) {
  //   // // Edge e1 = edgesVec.at(i)
  //   // Node * start = &vertices[findInArray(edges[i].start->id, vertices, numberOfNodes)];
  //   // Node * end = &vertices[findInArray(edges[i].end->id, vertices, numberOfNodes)];
  //   // // Edge * pred = findEdge(edges[i].pred->start->id,edges[i].pred->end->id,edgesVec);

  //   // Edge e1(end, edges[i].weight, start, edges[i].pred, edges[i].succ);
  //   // Edge e2(start, edges[i].weight, end, edges[i].pred, edges[i].succ);
  //   // uint32_t idx1 = findInArray(edges[i].start->id, vertices, numberOfNodes);
  //   // uint32_t idx2 = findInArray(edges[i].end->id, vertices, numberOfNodes);
  //   // this->adjacencyList[idx1].push_back(e1);
  //   // this->adjacencyList[idx2].push_back(e2);
  // }
}

Graph::~Graph() {
  delete[] adjacencyList;
  adjacencyList = nullptr;
  delete[] vertices;
  vertices = nullptr;
}

Graph dummyGraph() {
  return Graph(0,0,0);
}

void Graph::addEdge(uint32_t node1Id, uint32_t node2Id) {
  // We choose them, so that they're smaller than numberOfNodes
  const uint32_t weight = randomGen.generateRandomNumber(1, maxEdgeWeight);
  if(weight == 0)
    std::cout << "ERRRRRRRRRRROR\n";
  this->adjacencyList[node1Id].push_back(Edge(&vertices[node2Id], weight, &vertices[node1Id]));
  this->adjacencyList[node2Id].push_back(Edge(&vertices[node1Id], weight, &vertices[node2Id]));
}

void printAdajcencyList(std::vector<Edge>* adjList, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    std::cout << "node "<< i << std::endl;
    for (const Edge& e: adjList[i])
      std::cout << "edge to "<< e.end->id << " cost = " << e.weight << "; ";
    std::cout << std::endl;
  }
}

/*
Lists should have the same number of vectors.
@return returns true if they are different
*/
bool compareAdajcencyLists(std::vector<Edge>* adjList1, std::vector<Edge>* adjList2, uint32_t listSize) {
  for (uint32_t i = 0; i < listSize; ++i) {
    if (adjList1[i].size() != adjList2[i].size()) {
      std::cout << "Different number of elementsa at "<< i << std::endl;
      for (uint32_t j = 0; j < adjList1[i].size(); ++j)
        std::cout << adjList1[i].at(j).end->id << ",";
      std::cout << std::endl;
      for (uint32_t j = 0; j < adjList2[i].size(); ++j)
        std::cout << adjList2[i].at(j).end->id << ",";
      std::cout << std::endl;
      return true;
    } else {
      for (uint32_t j = 0; j < adjList1[i].size(); ++j) {
        if (adjList1[i].at(j).end != adjList2[i].at(j).end && adjList1[i].at(j).end != adjList2[i].at(j).end)
        {
          std::cout << "Difference in row " << i << " at element "<< j << " of value: "
            << adjList1[i].at(j).start->id << ","<< adjList1[i].at(j).end->id << " vs "
            << adjList2[i].at(j).start->id << ","<< adjList2[i].at(j).end->id  << std::endl;
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
      if(!edge.end->visited) {
        toVisit.push(edge.end);
        edge.end->visited = true;
      }
      //printVisitedStatus();
    }
  }
} 

void printEdge(Edge e) {
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

/*
copies adjlist of a graph to given copyOfAdjacencyList
*/
void Graph::copyAdjacencyList(std::vector<Edge>* copyOfAdjacencyList) {
  for (uint32_t i = 0; i < numberOfNodes; ++i) 
    for (const Edge& edge : adjacencyList[i]) {
      Node* start = &vertices[edge.start->id];  
      Node* end = &vertices[edge.end->id];      
      Edge copyOfEdge(end, edge.weight, start, edge.pred, edge.succ);
      copyOfAdjacencyList[i].push_back(copyOfEdge);
    }
}

/*
Resets modified priority_queue edge weights to those, present in adjacancyList of graph g
*/
// void resetQueue(std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparatorOnPointers>& toVisit, Graph * g) {
//   std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparatorOnPointers> tmpQueue;
//   std::swap(tmpQueue, toVisit);
//   while (!tmpQueue.empty()) { //unload queue to vector
//     Edge * edge = tmpQueue.top();
//     tmpQueue.pop();
//     Edge * e = g->findEdge(edge->start->id, edge->end->id);
//     //necessery for good behaviour of algorithm
//     // e->pred = edge->pred;
//     // e->succ = edge->succ;
//     edge->weight = e->weight;

//     //visited?

//     // std::cout << "Qunload found: \n";
//     // printEdge(*e);

//     toVisit.push(edge);
//     //TODO: should we check if it adds all edges correctly?
//   }
// }


//TODO: make sure generator, generates in range 1-something aka non zero
/*
Resets all non zero edges to edges from adjacencyList.
For Edges with weight == 0 we reset: visited(end and start), pred, succ
*/
void resetCopyOfAdjacencyList(std::vector<Edge>*& lacalCopyOfAdjacencyList, Graph * g) {
  for (uint32_t i = 0; i < g->numberOfNodes; ++i) {
    // std::cout << "lacalCopyOfAdjacencyList[" << i << "].size() " << lacalCopyOfAdjacencyList[i].size() << std::endl;
    // for (uint32_t jj = 0; jj < lacalCopyOfAdjacencyList[i].size(); ++jj) {
    //   std::cout << lacalCopyOfAdjacencyList[i].at(jj).end->id << ", ";
    // }
    // std::cout << std::endl;
    // std::cout << "g->adjacencyList[" << i << "].size() " << g->adjacencyList[i].size() << std::endl;
    // for (uint32_t jj = 0; jj < g->adjacencyList[i].size(); ++jj) {
    //   std::cout << g->adjacencyList[i].at(jj).end->id << ", ";
    // }
    // std::cout << std::endl;
    for (uint32_t j = 0; j < lacalCopyOfAdjacencyList[i].size(); ++j) {

      if (lacalCopyOfAdjacencyList[i].at(j).weight) { // if non zero edge - reset it to default (adjacencyList value)
        Edge * pred = lacalCopyOfAdjacencyList[i].at(j).pred;
        Edge * succ = lacalCopyOfAdjacencyList[i].at(j).succ;
        lacalCopyOfAdjacencyList[i].at(j).weight = g->adjacencyList[i].at(j).weight; //TODO: does it work as intended?
        lacalCopyOfAdjacencyList[i].at(j).pred = pred;
        lacalCopyOfAdjacencyList[i].at(j).succ = succ;
      } else { // weight == 0 aka edge in tree
        Edge * pred = &lacalCopyOfAdjacencyList[i].at(j);
        lacalCopyOfAdjacencyList[i].at(j).pred =pred;//nullptr // TODO: is this correct way to do this?
        lacalCopyOfAdjacencyList[i].at(j).succ = nullptr;
      }
      lacalCopyOfAdjacencyList[i].at(j).end->visited = false;
      lacalCopyOfAdjacencyList[i].at(j).start->visited = false;
    }
  }
}

/*
TODO: CHECK IF WORKS CORRECTLY
Return array of Edges with non empty fields passed to second parameter.
*/
void cropArraytoVector(Edge* * tmptreeEdges, std::vector<Edge>& vec, uint32_t treeSize) {
  for (uint32_t i = 0; i < treeSize; ++i)
    vec.push_back(*(tmptreeEdges)[i]);
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
Returns array of Nodes from array of Edges
*/
Node * getNodesFromEdges(std::vector<Edge> treeEdges, uint32_t treeSize) {
  std::cout << "getNodesFromEdges nodes: " << treeEdges.at(0).start->id <<std::endl;

  Node * nodeArray = new Node[treeSize + 1];

  std::cout << "HELLLO " <<std::endl;
  nodeArray[0] = *treeEdges.at(0).start;
  std::cout << "HELLLO " <<std::endl;
  // for (uint32_t i = 0; i < treeSize; ++i)
    // std::cout << "treeEdges[" << i << "] "<< treeEdges[i] <<std::endl;

  for (uint32_t i = 0; i < treeSize; ++i) {
    std::cout << "HELLLO " <<std::endl;
    std::cout << ", " << treeEdges.at(i).end->id <<std::endl;
    nodeArray[i + 1] = *treeEdges.at(i).end;
  }
  return nodeArray;
}

/*
Returns array of Nodes from array of Edges
*/
void getNodesFromEdges(Edge* treeEdges, Node** treeNodes, uint32_t treeSize) {
  std::cout << "getNodesFromEdges nodes: " << treeEdges[0].start->id;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  std::cout << "HELLLO " <<std::endl;
  (*treeNodes)[0] = *treeEdges[0].start;
  for (uint32_t i = 0; i < treeSize; ++i) {
    std::cout << ", " << treeEdges[i].end->id;
    (*treeNodes)[i + 1] = *treeEdges[i].end;
  }
}

/*
Finds both instances of edge and updates their start.
*/
void updateStart(Edge * e, Node * start, std::vector<Edge>* adjList) {
  e->start = start;
  findEdge(e->end->id, e->start->id, adjList)->start = start;
}
/*
Finds both instances of edge and updates their end.
*/
void updateEnd(Edge * e, Node * end, std::vector<Edge>* adjList) {
  e->end = end;
  findEdge(e->end->id, e->start->id, adjList)->end = end;
}
/*
Finds both instances of edge and updates their weights.
*/
void updateWeight(Edge * e, uint32_t weight, std::vector<Edge>* adjList) {
  e->weight = weight;
  findEdge(e->end->id, e->start->id, adjList)->weight = weight;
}

/*
Finds both instances of edge and updates their pred.
*/
void updatePred(Edge * e, Edge * pred, std::vector<Edge>* adjList) {
  e->pred = pred;
  findEdge(e->end->id, e->start->id, adjList)->pred = pred;
}

/*
Finds both instances of edge and creates pred loop.
*/
void updatePred(Edge * e, std::vector<Edge>* adjList) {
  e->pred = e;
  Edge * e2 = findEdge(e->end->id, e->start->id, adjList);
  e2->pred = e2;
}

/*
Finds both instances of edge and updates their succ.
*/
void updateSucc(Edge * e, Edge * succ, std::vector<Edge>* adjList) {
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

  std::vector<uint32_t> originalTerminals;
  for (uint32_t i = 0; i < terminals.size(); ++i)
    originalTerminals.push_back(terminals.at(i));

  std::cout << "dupa 1" << std::endl;
  std::cout << "terminalsWEJSCICE: ";
  for (uint32_t i = 0; i < terminals.size(); ++i)
    std::cout << terminals.at(i) << ", ";
  std::cout << std::endl;

  //init
  resetVisitedStatus();
  std::vector<Edge>* lacalCopyOfAdjacencyList = new std::vector<Edge>[numberOfNodes];
  copyAdjacencyList(lacalCopyOfAdjacencyList);
  Edge * tmptreeEdges[numberOfNodes - 1]; // no chyba nie tyle I guess ale maxymalnie tyle
      std::cout << "dupa 1" << std::endl;

  uint32_t treeSize = 0;
  bool foundTerminal = false;
  std::priority_queue<Edge *, std::vector<Edge *>, EdgeWeightComparatorOnPointers> toVisit;
      std::cout << "dupa 1" << std::endl;

  Edge dummyEdge(true);
  for (Edge& edge : lacalCopyOfAdjacencyList[terminals.at(0)]) {
    // updatePred(&edge, &dummyEdge, lacalCopyOfAdjacencyList);
    toVisit.push(&edge);
  }
  vertices[terminals.at(0)].visited = true;
  terminals.erase(terminals.begin());

  std::cout << "terminalsPOZA PIUERWSZEY: ";
  for (uint32_t i = 0; i < terminals.size(); ++i)
    std::cout << terminals.at(i) << ", ";
  std::cout << std::endl;



  do { //TMP NOTE - should be 9

    // std::cout << "adjacencyList" << std::endl;
    // printAdajcencyList(adjacencyList, numberOfNodes);
    // std::cout << "lacalCopyOfAdjacencyList" << std::endl;
    // printAdajcencyList(lacalCopyOfAdjacencyList, numberOfNodes);
    std::cout << "Comparing Adajcency Lists" << std::endl;
    if (compareAdajcencyLists(adjacencyList, lacalCopyOfAdjacencyList, numberOfNodes)) {
      std::cout << "Comparing Adajcency Lists Failed" << std::endl;
      return dummyGraph();
    }

    std::cout << "terminals "<< terminals.size() << std::endl;
    for (uint32_t i = 0; i < terminals.size(); ++i)
      std::cout << terminals.at(i) << ", ";
    std::cout << std::endl;




    if (foundTerminal) {
      std::cout << "dupaFIRTS 1" << std::endl;
      while(!toVisit.empty())
        toVisit.pop();
      std::cout << "dupaFIRTS 2" << std::endl;

      //TMP SOLUTION
      std::vector<Edge> tmpEdges;
      // cropArray(tmptreeEdges, &tmpEdges, treeSize);
      cropArraytoVector(tmptreeEdges, tmpEdges, treeSize);

      for (uint32_t i = 0; i < tmpEdges.size(); ++i)
        std::cout << "tmpEdges.at(" << i << ").end->id" << tmpEdges.at(i).end->id << "; tmpEdges.at(" << i << ").start->id" << tmpEdges.at(i).start->id << std::endl;

      std::cout << "dupaFIRTS 3" << std::endl;
      std::cout << "treeSize "<< treeSize << std::endl;

      // Node * tmpNodes = new Node[treeSize + 1];
      //       std::cout << "dupaFIRTS 4" << std::endl;

      Node * tmpNodes = getNodesFromEdges(tmpEdges, treeSize);
      std::cout << "dupaFIRTS 4 NODES" << std::endl;
      for (uint32_t i = 0; i < treeSize + 1; ++i) //for each node in current tree
      {
        std::cout << tmpNodes[i].id << ", ";
      }
      std::cout << std::endl;

      Edge* * e = new Edge*[treeSize + 1];
      for (uint32_t i = 0; i < treeSize + 1; ++i) //for each node in current tree
      {
        //TMP SOLUTION -find edges that are part of a tree
        std::cout << "dupaFIRTS 4.1 " << i << std::endl;
        for (Edge& edge : lacalCopyOfAdjacencyList[tmpNodes[i].id]) {
          if (!edge.weight)
          {
            e[i] = &edge;
            std::cout << "foudn e["<<i<<"].weight " << e[i]->weight << std::endl;

          }
        }

        std::cout << "dupaFIRTS 4.2" << std::endl;

        // add all other edges and set thier pred as edge from tree
        for (Edge& edge : lacalCopyOfAdjacencyList[tmpNodes[i].id]) {
          std::cout << "dupaFIRTS 4.2.0" << std::endl;

          if (edge.weight) // normal procedure
            updatePred(&edge, e[i], lacalCopyOfAdjacencyList);
          else // if it's 0 it's part of a tree
            updatePred(&edge, lacalCopyOfAdjacencyList);
          std::cout << "dupaFIRTS 4.2.1" << std::endl;
          std::cout << "edge.weight " << edge.weight << std::endl;
          std::cout << "ee["<<i<<"]->weight " <<e[i]->weight << std::endl;

          updateWeight(&edge, (edge.weight+e[i]->weight), lacalCopyOfAdjacencyList);
          toVisit.push(&edge);
          std::cout << "dupaFIRTS 4.2.2" << std::endl;

          std::cout << "adding:\n";
          printEdge(edge);
        }
        std::cout << "dupaFIRTS 4.3" << std::endl;

      }
      delete[] tmpNodes;
      tmpNodes = nullptr;
      delete[] e;
      e = nullptr;
      std::cout << "dupaFIRTS 5" << std::endl;
      foundTerminal = false;
    }
    std::cout << "dupa3.1" << std::endl;





    //main loop
    while(!foundTerminal) { // O(n^2)
          std::cout << "TERMINALS "<< terminals.size() << std::endl;

      //remoeve all edges that lead to already visited nodes
      if (!toVisit.empty())
        while (toVisit.top()->end->visited)
          toVisit.pop();

      std::cout << "dupa3.1.5" << std::endl;

      if (toVisit.empty()) {
        std::cout << "End in loop Dijkstra; dummy graph" << std::endl;
        return dummyGraph();
      }
      std::cout << "dupa3.2" << std::endl;

      Edge * e = toVisit.top();
      toVisit.pop();
      e->end->visited = true;
      uint32_t nextNodeIndex = e->end->id;
      std::cout << "dupa3.3" << std::endl;



      // check if we have found terminal
      int32_t idx = findInVector(nextNodeIndex, terminals);
      if (idx > -1) {
        foundTerminal = true;
        std::cout << "dupa3.4" << std::endl;

        //remove found terminal from the list
        terminals.erase(terminals.begin() + idx);
        std::cout << "dupa3.5" << std::endl;

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
        std::cout << "dupa3.6" << std::endl;
        if (e->pred == nullptr) {
          std::cerr << "WARNING, e->pred is null poiter: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
          tmptreeEdges[treeSize] = findEdge(e->start->id, e->end->id, adjacencyList);
          std::cout << "dupa3.6.01" << std::endl;

          if (tmptreeEdges[treeSize] == nullptr) {
            std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
            return dummyGraph();
          }
          std::cout << "dupa3.6.02" << std::endl;
          printEdge(*e);

          ++treeSize;
          updateWeight(e, 0, lacalCopyOfAdjacencyList);
          updatePred(e, lacalCopyOfAdjacencyList);

          // return dummyGraph();
        }
        //TODO: clean code below?
        // printEdge(*e);
        // zero edges and add their original instance from adjacecnyList to tmptreeEdgeas
        else {
          std::cerr << "e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
          std::cerr << "e->pred->start->id: "<< e->pred->start->id << "; e->pred->end->id: " << e->pred->end->id << "; e->weight: "<< e->pred->weight << std::endl;
        }
        while (e->pred != nullptr && e->pred->start->id != e->start->id) {
          std::cout << "dupa3.6.1" << std::endl;
          tmptreeEdges[treeSize] = findEdge(e->start->id, e->end->id, adjacencyList);
          std::cout << "dupa3.6.2" << std::endl;
          if (tmptreeEdges[treeSize] == nullptr) {
            std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
            return dummyGraph();
          }
          std::cout << "dupa3.6.3" << std::endl;
          printEdge(*e);

          ++treeSize;
          updateWeight(e, 0, lacalCopyOfAdjacencyList);
          Edge * oldPred = e->pred;
          updatePred(e, lacalCopyOfAdjacencyList);

          std::cout << "dupa3.6.4" << std::endl;

          std::cerr << "e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
          if (oldPred == nullptr) {
            std::cerr << "huj" << std::endl;
            break;
          }
          else {
            std::cerr << "oldPred->start->id: "<< oldPred->start->id << "; oldPred->end->id: " << oldPred->end->id << "; e->weight: "<< oldPred->weight << std::endl;
            e = oldPred;
          }
        } // untill we reach beginning of the path
        std::cout << "dupa3.7" << std::endl;
          printEdge(*e);

        //add last edge
        // tmptreeEdges[treeSize] = findEdge(e->start->id, e->end->id, adjacencyList);
        // if (tmptreeEdges[treeSize] == nullptr) {
        //   std::cerr << "error, edge not found: e->start->id: "<< e->start->id << "; e->end->id: " << e->end->id << "; e->weight: "<< e->weight << std::endl;
        //   return dummyGraph();
        // }
        // ++treeSize;
        updateWeight(e, 0, lacalCopyOfAdjacencyList);

        std::cout << "dupa3.8" << std::endl;

        resetCopyOfAdjacencyList(lacalCopyOfAdjacencyList, this);
        std::cout << "dupa3.9" << std::endl;
        break;
      }
      for (Edge& edge : lacalCopyOfAdjacencyList[nextNodeIndex]) {
        if (!edge.end->visited) {
          if (edge.weight) // normal procedure
            // edge.pred = e;
            updatePred(&edge, e, lacalCopyOfAdjacencyList);
          else // if it's 0 it's part of a tree
            updatePred(&edge, lacalCopyOfAdjacencyList);
            // edge.pred = &edge;
          updateWeight(&edge, (edge.weight + e->weight), lacalCopyOfAdjacencyList);
          // edge.weight += e->weight;
          toVisit.push(&edge);
          std::cout << "adding:\n";
          printEdge(edge);
        }
      }
      std::cout << "dupa4" << std::endl;
      // printData();
      for (uint32_t i = 0; i < treeSize; ++i)
        std::cout << tmptreeEdges[i]->start->id << "->" <<  tmptreeEdges[i]->end->id << "; ";
      std::cout << std::endl;
      std::cout << "treeSize: " << treeSize << std::endl;
      std::cout << "terminals: ";
      for (uint32_t i = 0; i < terminals.size(); ++i)
        std::cout << terminals.at(i) << ", ";
      std::cout << std::endl;
    }
  } while(!terminals.empty()); // O(k)

  //ISSUES, MULTIPLE SAME EDGES, DOESNT SEE MORE TERMINALS

  //TODO: print edges in separate function
  // print edges
  if (printFlag && numberOfNodes < 1000) {
    for (uint32_t i = 0; i < treeSize; ++i)
      std::cout << tmptreeEdges[i]->start->id << "->" <<  tmptreeEdges[i]->end->id << "; ";
    std::cout << std::endl;
  }
  std::cout << "Original terminals: ";
  for (uint32_t i = 0; i < originalTerminals.size(); ++i)
    std::cout << originalTerminals.at(i) << ", ";
  uint64_t totalWeight = 0;
  for (uint32_t i = 0; i < treeSize; ++i)
    totalWeight += tmptreeEdges[i]->weight;
  std::cout << "Dijkstra totalWeight = " << totalWeight << std::endl;

  std::vector<Edge> treeEdgesVec;
  cropArraytoVector(tmptreeEdges, treeEdgesVec, treeSize);

  Edge * treeEdges = new Edge[treeSize];
  for (uint32_t i = 0; i < treeEdgesVec.size(); ++i) {
    treeEdges[i] = treeEdgesVec.at(i);
  }
  // cropArray(tmptreeEdges, &treeEdges, treeSize);
  Node * treeNodes = getNodesFromEdges(treeEdgesVec, treeSize);
    std::cout << "dupa 123" << std::endl;

  // getNodesFromEdges(treeEdges, &treeNodes, treeSize);
  Graph steinerTree = Graph(treeEdges, treeNodes, static_cast<uint32_t>(treeSize + 1), treeSize, printFlag);
  std::cout << "dupa 321" << std::endl;

  delete[] lacalCopyOfAdjacencyList;
  lacalCopyOfAdjacencyList = nullptr;
  std::cout << "dupa 321" << std::endl;

  delete[] treeEdges;
  treeEdges = nullptr;
  std::cout << "dupa 321" << std::endl;

  delete[] treeNodes;
  treeNodes = nullptr;
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
  for (Edge& edge : adjacencyList[0]) //init 
    toVisit.push(&edge);
  vertices[0].visited = true;
  
  Edge * tmptreeEdges[numberOfNodes - 1];
  uint64_t totalWeight = 0;
  for (uint32_t treeSize = 0; treeSize < numberOfNodes - 1; ++treeSize)
  {
    while (toVisit.top()->end->visited) //remoeve all edges that lead to already visited nodes
      toVisit.pop();

    if (toVisit.size() == 0) 
    {
      std::cout << "End in loop PimMST; dummy graph" << std::endl;
      return dummyGraph();
    }
    Edge * e = toVisit.top();
    tmptreeEdges[treeSize] = e;
    toVisit.pop();
    e->end->visited = true;
    uint32_t nextNodeIndex = e->end->id;
    for (Edge& edge : adjacencyList[nextNodeIndex])
    {
      if (!edge.end->visited)
      {
        toVisit.push(&edge);
      }
    }
    totalWeight += e->weight;
  }
  if (printFlag && numberOfNodes < 1000)
  {
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
    for (const auto& edge : adjacencyList[i]) {
      std::cout << static_cast<uint32_t>(edge.end->id) << " cost(" << edge.weight << ") "; 
      if(edge.weight == 0)
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
  for (const auto& edge : adjacencyList[node1Id])
    if (edge.end->id == node2Id)
      return true;
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
}



#ifndef HELPERS_HPP
#define HELPERS_HPP

/*  ADJACENCY LIST RELATED  */
/*
  Compares 2 adjacency lists and prints where they differ if they differ.
  Returns 0 if all is fine.
*/
bool compareAdajcencyLists(std::vector<std::shared_ptr<Edge>>* adjList1, std::vector<std::shared_ptr<Edge>>* adjList2, uint32_t listSize);
/*
  Given a graph it creates a copy of it's adjacency list and saves it in secend function argument.
  Edges are ceated anew, there is no passing of Edge pointers.
*/
void copyAdjacencyListFromGraphWithNewNodeInstances(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList);
/*
  Moves all nonempty edges to the front of given adjacency list.
*/
void sortAdjacencyList(uint32_t currentAdjListSize, std::vector<std::shared_ptr<Edge>>*& adjList);
/*
  Converts adjacency list into a vector of edges with no repeats
  e.g. if we have an edge 1-2 then edge 2-1 won't be added
*/
std::vector<std::shared_ptr<Edge>> flattenAdjacencyList(
    std::vector<std::shared_ptr<Edge>>* adjacencyList,
    uint32_t numberOfNodes);
/*  RESETERS  */
/*
  Resets visited status and weight in copy of adjacency list.
*/
void resetVisitedStatusAndWeightInCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList, std::shared_ptr<Graph> g);
/*
  Resets all edge fields in copy of adjacency list.
*/
void fullResetCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList, std::shared_ptr<Graph> g);

/*  FINDERS  */
/*
  Finds value in uint32_t vector.
  Returs either it's index inside this vector or -1 if value was not found.
*/
int32_t findInUintVector(uint32_t value, std::vector<uint32_t> vec);
/*
  Finds edge of given start and end in edge vector.
  Returs either it's index inside this vector or -1 if value was not found.
*/
int32_t findInEdgeVector(uint32_t start, uint32_t end, std::vector<std::shared_ptr<Edge>> vec);
/*
  Finds edge of given start and end in edge vector.
  Returs pointer to that edge or nullptr if edge was not found.
*/
std::shared_ptr<Edge> findInEdgeVectorAndReturnValue(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>> vec);
/*
  Finds value in uint32_t array.
  Returs either it's index or -1 if value was not found.
*/
int32_t findInArray(uint32_t value, std::shared_ptr<Node>* array, uint32_t arraySize);
/*
  Finds PseudoEdge of given start and end in PseudoEdge vector.
  Returs either it's index inside this vector or -1 if PseudoEdge was not found.
*/
int32_t findInTmpTree(uint32_t edgeStart, uint32_t edgeEnd, std::vector<PseudoEdge> vec);
/*
  Finds edge of given start and end in given adjacency list.
  Returs pointer to that edge or nullptr if edge w not found.
*/
std::shared_ptr<Edge> findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<std::shared_ptr<Edge>>* adjacencyList);
/*
  Returns instance of edge of weight equal to 0 adjacenct to node of given index.
*/
std::shared_ptr<Edge> findZeroEdgeInAdjacentTo(uint32_t idx, std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList);
/*
  Given a vector of shortest paths it returns weight of the one going from startId to endId (or from endId to startId).
*/
uint32_t findShortestPathAndReturnWeight(
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<Edge>>>> shortestPathEdgeVec,
    uint32_t startId,
    uint32_t endId);

/*  UPDATERS  */
/*
  Updates weight of given edge and it's mirror image (updates e.g. 1->2 and 2->1).
*/
void updateWeight(std::shared_ptr<Edge>& e, uint32_t weight, std::vector<std::shared_ptr<Edge>>*& adjList);
/*
  Updates predecessor of given edge and it's mirror image (updates e.g. 1->2 and 2->1).
*/
void updatePred(std::shared_ptr<Edge>& e, std::shared_ptr<Edge> pred, std::vector<std::shared_ptr<Edge>>* adjList);
/*
  Updates predecessor of given edge and it's mirror image (updates e.g. 1->2 and 2->1) to self loop.
*/
void updatePredToLoop(std::shared_ptr<Edge>& e, std::vector<std::shared_ptr<Edge>>* adjList);

/*  PSEUDO-EDGE RELATED  */
/*
  Returns number of unique numbers in given PseudoEdge vector.
*/
uint32_t tmpPseudoEdgeFindUniqueNumbers(std::vector<PseudoEdge> vec);
/*
  Returns unique numbers in given PseudoEdge vector.
*/
std::vector<uint32_t> tmpPseudoEdgeReturnUniqueNumbers(std::vector<PseudoEdge> vec);
/*
  Adds edge to PseudoEdge vector.
*/
void addTotmptreeEdgesIfNotAlreadyIn(std::vector<PseudoEdge>& tmptreeEdges, std::shared_ptr<Edge> edg);

/*  ADDERS  */
/*
  Adds node to given vector if it is not already in and increases counter.
*/
void addNodeIfNotAlreadyIn(uint32_t nodeId, std::vector<uint32_t>& vec, uint32_t& counter);
/*
  Adds edge to given vector if it is not already in and increases counter.
*/
void addEdgeIfNotAlreadyIn(std::shared_ptr<Edge> edge, std::vector<std::shared_ptr<Edge>>& vec, uint32_t& counter);
/*
  Updates weight of given edge to 0.
  Updates its predecessor to loop.
  Finds this edge in given original adjacency list.
  Adds it to a tmpTree of PseudoEdges.
*/
void embedEdgeIntoTree(
    std::shared_ptr<Edge> &e,
    std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
    std::vector<std::shared_ptr<Edge>>* &adjacencyList,
    std::vector<PseudoEdge> &tmptreeEdges);

/*  OTHER  */
/*
  Given array of type Node returns vector of type uint32_t.
*/
std::shared_ptr<std::vector<uint32_t>> nodeArrayToUint(std::shared_ptr<Node>* vertices, uint32_t numberOfNodes);
/*
  Returns number of edges in Complete graph.
*/
uint32_t numberOfEdgesInClique(uint32_t numberOfNodes);
/*
  Returns weight of given path.
*/
uint32_t calculatePathWeight(std::vector<std::shared_ptr<Edge>> path);
/*
  Checks if nth bit of given number is equal to 1.
*/
bool isNthBitSet(uint32_t num, uint32_t n);

#endif

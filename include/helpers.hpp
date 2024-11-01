#ifndef HELPERS_HPP
#define HELPERS_HPP


bool compareAdajcencyLists(std::vector<std::shared_ptr<Edge>>* adjList1, std::vector<std::shared_ptr<Edge>>* adjList2, uint32_t listSize);
// void copyAdjacencyListFromGraph(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList);
void copyAdjacencyListFromGraphWithNewNodeInstances(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<Edge>>*& copyOfAdjacencyList);
int32_t findInUintVector(uint32_t value, std::vector<uint32_t> vec);
int32_t findInEdgeVector(uint32_t start, uint32_t end, std::vector<std::shared_ptr<Edge>> vec);
std::shared_ptr<Edge> findInEdgeVectorAndReturnValue(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>> vec);
int32_t findInArray(uint32_t value, std::shared_ptr<Node>* array, uint32_t arraySize);
int32_t findInTmpTree(uint32_t edgeStart, uint32_t edgeEnd, std::vector<PseudoEdge> vec);
std::shared_ptr<Edge> findEdge(uint32_t edgeStart, uint32_t edgeEnd, std::vector<std::shared_ptr<Edge>>* adjacencyList);
void resetCopyOfAdjacencyList(std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList, std::shared_ptr<Graph> g);
void updateWeight(std::shared_ptr<Edge>& e, uint32_t weight, std::vector<std::shared_ptr<Edge>>*& adjList);
void updatePred(std::shared_ptr<Edge>& e, std::shared_ptr<Edge> pred, std::vector<std::shared_ptr<Edge>>* adjList);
void updatePredToLoop(std::shared_ptr<Edge>& e, std::vector<std::shared_ptr<Edge>>* adjList);
uint32_t tmpPseudoEdgeFindUniqueNumbers(std::vector<PseudoEdge> vec);
std::vector<uint32_t> tmpPseudoEdgeReturnUniqueNumbers(std::vector<PseudoEdge> vec);
void addTotmptreeEdgesIfNotAlreadyIn(std::vector<PseudoEdge>& tmptreeEdges, std::shared_ptr<Edge> edg);
std::shared_ptr<Edge> findZeroEdgeInAdjacentTo(uint32_t idx, std::vector<std::shared_ptr<Edge>>*& localCopyOfAdjacencyList);
uint32_t numberOfEdgesInClique(uint32_t numberOfNodes);
void removeSingleEdgeFromAdjacencyList(std::shared_ptr<Node> start, std::shared_ptr<Node> end, std::vector<std::shared_ptr<Edge>>*& adjList);
void removeEdgeFromAdjacencyList(std::shared_ptr<Edge> e, std::vector<std::shared_ptr<Edge>>*& adjList);
void sortAdjacencyList(uint32_t currentAdjListSize, std::vector<std::shared_ptr<Edge>>*& adjList);
void addNodeIfNotAlreadyIn(uint32_t nodeId, std::vector<uint32_t>& vec, uint32_t& counter);
void addEdgeIfNotAlreadyIn(std::shared_ptr<Edge> edge, std::vector<std::shared_ptr<Edge>>& vec, uint32_t& counter);
void embedEdgeIntoTree(
    std::shared_ptr<Edge> &e,
    std::vector<std::shared_ptr<Edge>>* &localCopyOfAdjacencyList,
    std::vector<std::shared_ptr<Edge>>* &adjacencyList,
    std::vector<PseudoEdge> &tmptreeEdges);

#endif

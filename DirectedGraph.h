#ifndef DirectedGraph_h
#define DirectedGraph_h

#include <vector>
#include <iostream>

using matrix2d = std::vector<std::vector<int>> ;

/// Class representing directed graph. By default graph like this one does not
/// have weights, but they can be generated. This class provides methods that
/// change graph representation and perform some simple algorithms. 
class DirectedGraph
{
public:
	/// Public contructors ///

	/// Most basic directed graph constructor. Its first argument represents 2d matrix that
	/// stores information about graph in one of three possible representations
	/// Second argument contains information about representation in a form of string. Accepted 
	/// input : 
	/// - "adjacency_list" - graph is stored in a form of adjacency list
	/// - "adjacency_matrix" - graph is stored in a form of adjacency matrix
	/// - "incidence_matrix" - graph is stored as incidence matrix
	/// Created graph has no weights, so its distance and p matrices are set to default values
	DirectedGraph (matrix2d data, std::string representation);

	/// This static metod returns random graph generated using information from arguments. 
	/// First arguments indicates how many vertices generated graph will have, and second one
	/// gives probability of edge existence between every pair of vertices. 
	/// So for example if we have 3 vertices, there are 6 posible edges - each can be created
	/// with given probability.
	static DirectedGraph* GenerateRandomGraphBasedOnProbability(unsigned numberOfVertices, double probability);

	/// This function generates random strongly connected graph. This graph can have anything between
	/// 1 and 10 vertices - more would make plot hard to read. 
	static DirectedGraph* GenerateRandomStronglyConnectedGraph();

	///


	/// Getters ///

	/// Calculates and returns number of edges in a graph
	unsigned GetNumberOfEdges() const;

	/// Returns number of vertices in a graphs
	unsigned GetNumberOfVertices() const { return m_numberOfVertices; }

	/// Returns boolean that represents information if 
	/// weights have been generated for stored graph
	bool HasWeights() const { return m_hasWeights; }

	///


	/// Representation changing ///

	/// Changes graph representation to adjacency matrix
	/// If graph is already represented by adjacency matrix, 
	/// then it checks it and quits
	void ChangeToAdjacencyMatrix();

	/// Changes graph representation to adjacency list.
	void ChangeToAdjacencyList();

	/// Changes graph representation to incidence matrix
	void ChangeToIncidenceMatrix();

	///


	/// Printing and drawing ///

	/// Prints graph in currently stored representation
	void PrintGraph(std::ostream& placeToPrint) const;

	/// Prints graphs weight matrix if one has been generated. 
	/// Otherwise, it prints suitable information and returns
	void PrintWeightMatrix(std::ostream& placeToPrint) const;

	///


	/// Algorithms ///

	/// Generates graph weights form range [lowerBound, upperBound]
	/// Every graph edge has a weight after this algorithm is used.
	/// Sets m_hasWeights flag to true and allows to print weight
	/// matrix. 
	void GenerateWeights(int lowerBound, int upperBound);

	/// Implemetation of Kosaraju algorithm
	std::vector<int> Kosaraju();

	/// Implementation of Bellman - Ford algorithm
	bool BellmanFord(int source, bool print = true);

	/// Implementation of Johnsons algorithm
	matrix2d Johnson();

	/// Implementation of Pagerank a) version
	void PageRankA();

	/// Implementation of Pagerank b) version
	void PageRankB();

	///
	

private: /// fields

	enum class GraphRepresentation
	{
		AdjacencyList = 1,
		AdjacencyMatrix = 2,
		IncidenceMatrix = 3
	};

	GraphRepresentation m_currentRepresentation;

	bool m_isStronglyConected;

	bool m_hasWeights;

	unsigned m_numberOfVertices;

	matrix2d m_encodedGraphData;	

	matrix2d m_weightMatrix;

	matrix2d distances;

	matrix2d p;


private: /// methods

	/// Representation changing ///

	void AdjacencyListToAdjacencyMatrix();

	void IncidenceMatrixToAdjacencyMatrix();

	void AdjacencyMatrixToAdjacencyList();

	void IncidenceMatrixToAdjacencyList();

	void AdjacencyMatrixToIncidenceMatrix();

	void AdjacencyListToIncidenceMatrix();


	/// Algorithms ///

	void DFSvisit(int vertexID, std::vector<int>& visitTime, std::vector<int>& processTime, int& time);

	void ComponentsR(int nr, int vertexID, matrix2d transposedGraph, std::vector<int>& comp);

	DirectedGraph* Add_S();

	void Init(int source);

	void Relax(matrix2d& weights, int u, int v, int source);

	void Dijkstra(matrix2d& weights, int source);

	void InitializeDistancesAndP();


	/// Private constructor ///

	DirectedGraph() : m_numberOfVertices(0), m_hasWeights(false), m_isStronglyConected(false) {}

};

#endif

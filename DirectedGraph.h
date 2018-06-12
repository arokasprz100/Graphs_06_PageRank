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

	/// Calculates and returns number of edges in a graph
	unsigned GetNumberOfEdges() const;

	/// Returns number of vertices in a graphs
	unsigned GetNumberOfVertices() const { return m_numberOfVertices; }

	/// Changes graph representation to adjacency matrix
	/// If graph is already represented by adjacency matrix, 
	/// then it checks it and quits
	void ChangeToAdjacencyMatrix();

	/// Changes graph representation to adjacency list.
	void ChangeToAdjacencyList();

	/// Changes graph representation to incidence matrix
	void ChangeToIncidenceMatrix();

	/// Prints graph in currently stored representation
	void PrintGraph(std::ostream& placeToPrint) const;

	/// Implementation of Pagerank a) version
	void PageRankA();

	/// Implementation of Pagerank b) version
	void PageRankB();	

private: /// fields

	enum class GraphRepresentation
	{
		AdjacencyList = 1,
		AdjacencyMatrix = 2,
		IncidenceMatrix = 3
	};

	GraphRepresentation m_currentRepresentation;

	bool m_isStronglyConected;

	unsigned m_numberOfVertices;

	matrix2d m_encodedGraphData;	

private: /// methods

	void AdjacencyListToAdjacencyMatrix();

	void IncidenceMatrixToAdjacencyMatrix();

	void AdjacencyMatrixToAdjacencyList();

	void IncidenceMatrixToAdjacencyList();

	void AdjacencyMatrixToIncidenceMatrix();

	void AdjacencyListToIncidenceMatrix();

	void ComponentsR(int nr, int vertexID, matrix2d transposedGraph, std::vector<int>& comp);

	void DFSvisit(int vertexID, std::vector<int>& visitTime, std::vector<int>& processTime, int& time);

	std::vector<int> Kosaraju();

	DirectedGraph() : m_numberOfVertices(0), m_isStronglyConected(false) {}

};

#endif

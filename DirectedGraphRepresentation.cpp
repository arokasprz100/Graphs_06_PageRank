#include "DirectedGraph.h"

void DirectedGraph::ChangeToAdjacencyMatrix()
{
	if (m_currentRepresentation == GraphRepresentation::AdjacencyMatrix)
		return;
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyList)
		AdjacencyListToAdjacencyMatrix();
	else if (m_currentRepresentation == GraphRepresentation::IncidenceMatrix)
		IncidenceMatrixToAdjacencyMatrix();

}

void DirectedGraph::AdjacencyListToAdjacencyMatrix()
{
	matrix2d newData(m_encodedGraphData.size(), std::vector<int>(m_encodedGraphData.size(), 0));

	for (unsigned i = 0; i<m_numberOfVertices; ++i)
		for (unsigned j = 0; j<m_encodedGraphData.at(i).size(); ++j)
			newData.at(i).at(m_encodedGraphData.at(i).at(j)) = 1;

	m_encodedGraphData = newData;
	m_currentRepresentation = GraphRepresentation::AdjacencyMatrix;
}

void DirectedGraph::IncidenceMatrixToAdjacencyMatrix()
{
	matrix2d newData(m_encodedGraphData.size(), std::vector<int>(m_encodedGraphData.size(), 0));
	unsigned edgeStart = 0;
	unsigned edgeEnd = 0;

	for (unsigned i = 0; i<m_encodedGraphData.at(0).size(); ++i)
	{
		for (unsigned j = 0; j<m_encodedGraphData.size(); ++j)
		{

			if (m_encodedGraphData.at(j).at(i) == -1)
				edgeStart = j;
			else if (m_encodedGraphData.at(j).at(i) == 1)
				edgeEnd = j;
		}
		newData.at(edgeStart).at(edgeEnd) = 1;
	}

	m_encodedGraphData = newData;
	m_currentRepresentation = GraphRepresentation::AdjacencyMatrix;
}

void DirectedGraph::ChangeToAdjacencyList()
{
	if (m_currentRepresentation == GraphRepresentation::AdjacencyList)
		return;
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyMatrix)
		AdjacencyMatrixToAdjacencyList();
	else if (m_currentRepresentation == GraphRepresentation::IncidenceMatrix)
		IncidenceMatrixToAdjacencyList();

}

void DirectedGraph::AdjacencyMatrixToAdjacencyList()
{
	matrix2d newData;
	for (unsigned i = 0; i<m_numberOfVertices; ++i)
	{
		newData.emplace_back();
		for (unsigned j = 0; j<m_numberOfVertices; ++j)
			if (m_encodedGraphData.at(i).at(j) == 1)
				newData.at(i).push_back(j);
	}

	m_encodedGraphData = newData;
	m_currentRepresentation = GraphRepresentation::AdjacencyList;
}

void DirectedGraph::IncidenceMatrixToAdjacencyList()
{
	matrix2d newData(m_encodedGraphData.size(), std::vector<int>());
	unsigned edgeStart = 0;
	unsigned edgeEnd = 0;

	for (unsigned i = 0; i<m_encodedGraphData.at(0).size(); ++i)
	{
		for (unsigned j = 0; j<m_encodedGraphData.size(); ++j)
		{
			if (m_encodedGraphData.at(j).at(i) == -1)
				edgeStart = j;
			else if (m_encodedGraphData.at(j).at(i) == 1)
				edgeEnd = j;
		}
		newData.at(edgeStart).push_back(edgeEnd);
	}

	m_encodedGraphData = newData;
	m_currentRepresentation = GraphRepresentation::AdjacencyList;

}

void DirectedGraph::ChangeToIncidenceMatrix()
{
	if (m_currentRepresentation == GraphRepresentation::IncidenceMatrix)
		return;
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyMatrix)
		AdjacencyMatrixToIncidenceMatrix();
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyList)
		AdjacencyListToIncidenceMatrix();
}

void DirectedGraph::AdjacencyMatrixToIncidenceMatrix()
{
	matrix2d incidenceMatrix(m_encodedGraphData.size(), std::vector<int>(GetNumberOfEdges(), 0));
	int edgeNumber = 0;

	for (unsigned i = 0; i < m_numberOfVertices; ++i)
	{
		for (unsigned j = 0; j < m_encodedGraphData.at(i).size(); ++j)
		{
			if (m_encodedGraphData[i][j] == 1)
			{
				incidenceMatrix[i][edgeNumber] = -1;
				incidenceMatrix[j][edgeNumber] = 1;
				edgeNumber++;
			}
		}
	}

	m_encodedGraphData = incidenceMatrix;
	m_currentRepresentation = GraphRepresentation::IncidenceMatrix;
}

void DirectedGraph::AdjacencyListToIncidenceMatrix()
{
	AdjacencyListToAdjacencyMatrix();
	AdjacencyMatrixToIncidenceMatrix();
}
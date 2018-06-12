#include "DirectedGraph.h"
#include <random>
#include <exception>
#include <set>

DirectedGraph::DirectedGraph (matrix2d data, std::string representation)
{
	m_encodedGraphData = data;
	if (representation == "adjacency_list")
		m_currentRepresentation = GraphRepresentation::AdjacencyList;
	else if (representation == "adjacency_matrix")
		m_currentRepresentation = GraphRepresentation::AdjacencyMatrix;
	else if (representation == "incidence_matrix")
		m_currentRepresentation = GraphRepresentation::IncidenceMatrix;

	m_numberOfVertices = data.size();
	m_isStronglyConected = false;

}

DirectedGraph* DirectedGraph::GenerateRandomGraphBasedOnProbability(unsigned numberOfVertices, double probability)
{
	DirectedGraph* graph = new DirectedGraph;
	graph->m_currentRepresentation = GraphRepresentation::AdjacencyMatrix;
	graph->m_numberOfVertices = numberOfVertices;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 100.0);

	for (unsigned i = 0; i<numberOfVertices; ++i)
	{
		graph->m_encodedGraphData.emplace_back();
		for (unsigned j = 0; j<numberOfVertices; ++j)
		{
			if (i == j)
				graph->m_encodedGraphData.at(i).push_back(0);
			else if (dis(gen) < probability)
				graph->m_encodedGraphData.at(i).push_back(1);
			else
				graph->m_encodedGraphData.at(i).push_back(0);
		}
	}

	return graph;
}

DirectedGraph* DirectedGraph::GenerateRandomStronglyConnectedGraph()
{
	std::random_device randDev;
	std::mt19937 gen(randDev());
	std::uniform_int_distribution<> ver(2, 10);
	std::uniform_int_distribution<> prob(1, 100);

	DirectedGraph * graph = DirectedGraph::GenerateRandomGraphBasedOnProbability(2, 1);

	while (graph->m_isStronglyConected == false)
	{
		delete graph;
		graph = DirectedGraph::GenerateRandomGraphBasedOnProbability(ver(gen), prob(gen));
		graph->Kosaraju();
	}
	return graph;
}


void DirectedGraph::PrintGraph (std::ostream& placeToPrint) const
{
	placeToPrint << std::endl;
	if (m_currentRepresentation == GraphRepresentation::AdjacencyList)
		placeToPrint<<"Adjacency List"<<std::endl<<std::endl;
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyMatrix)
		placeToPrint<<"Adjacency Matrix"<<std::endl<<std::endl;
	else if (m_currentRepresentation == GraphRepresentation::IncidenceMatrix)
		placeToPrint<<"Incidence Matrix"<<std::endl<<std::endl;

	for (unsigned i = 0; i<m_numberOfVertices; ++i)
	{
		placeToPrint << i << ": ";
		for (unsigned j =0; j<m_encodedGraphData.at(i).size(); ++j)
			placeToPrint<<m_encodedGraphData.at(i).at(j)<<" ";
		placeToPrint<<std::endl;
	}
	placeToPrint << std::endl;
}


unsigned DirectedGraph::GetNumberOfEdges() const
{
	unsigned counter = 0;
	if (m_currentRepresentation == GraphRepresentation::AdjacencyList)
	{
		
		for (unsigned i = 0; i<m_numberOfVertices; ++i)
			counter += m_encodedGraphData.at(i).size();
		return counter;
	}
	else if (m_currentRepresentation == GraphRepresentation::AdjacencyMatrix )
	{
		for (unsigned i =0; i<m_numberOfVertices; ++i)
			for (unsigned j =0; j<m_numberOfVertices; ++j)
				if (m_encodedGraphData.at(i).at(j) == 1)
					counter++;
		return counter;
	}
	else 
		return m_encodedGraphData.at(0).size();
}




std::vector<int> DirectedGraph::Kosaraju()
{
	ChangeToAdjacencyMatrix();
	std::vector<int> d(m_encodedGraphData.size(), -1);
	std::vector<int> f(m_encodedGraphData.size(), -1);
	int t=0;
	for(unsigned i=0; i<m_encodedGraphData.size(); ++i)
	{
		if(d[i] == -1)
			DFSvisit(i, d, f, t);
	}

	matrix2d transposedGraph (m_encodedGraphData.size(), std::vector<int>(m_encodedGraphData.size(), 0));

	for (unsigned i = 0; i < m_encodedGraphData.size(); ++i)
	{
		for (unsigned j = 0; j < m_encodedGraphData.at(i).size(); ++j)
		{
			if(m_encodedGraphData.at(i).at(j) == 1)
				transposedGraph[j][i]=1;
		}
	}

	int nr=0;
	std::vector<int> comp(m_encodedGraphData.size(), -1);
	std::vector<int> decreaseF(m_encodedGraphData.size(), 0);

	for(unsigned i = 0; i < m_encodedGraphData.size(); ++i)
		decreaseF[i]=i;

	for(unsigned i = 0; i < (m_encodedGraphData.size()-1); ++i)
	{
		for(unsigned j = (m_encodedGraphData.size()-1); j >= 1+i; --j)
		{
			if(f[j] > f[j-1])
			{
				int temp_f = f[j];
				f[j] = f[j-1];
				f[j-1] = temp_f;
				int temp = decreaseF[j];
				decreaseF[j]=decreaseF[j-1];
				decreaseF[j-1]=temp;
			}
		}
	}

	for (unsigned i = 0; i < comp.size(); ++i)
	{
		if(comp[decreaseF[i]] == -1)
		{
			nr=nr+1;
			comp[decreaseF[i]]=nr;
			ComponentsR(nr, decreaseF[i], transposedGraph, comp);
		}
	}
	
	 if(nr == 1)
		m_isStronglyConected = true;
	
	return comp;
}


void DirectedGraph::DFSvisit(int vertexID, std::vector<int>& visitTime, std::vector<int>& processTime, int& time)
{
	time=time+1;
	visitTime[vertexID]=time;
	for(unsigned i=0; i<m_encodedGraphData.size(); ++i)
	{
		if(m_encodedGraphData.at(vertexID).at(i) == 1)
		{
			if(visitTime[i] == -1)
				DFSvisit(i, visitTime, processTime, time);
		}
	}
	time=time+1;
	processTime[vertexID]=time;
}


void DirectedGraph::ComponentsR(int nr, int vertexID, matrix2d transposedGraph, std::vector<int>& comp)
{
	for(unsigned i=0; i<transposedGraph.size(); ++i)
	{
		if(transposedGraph.at(vertexID).at(i) == 1)
		{
			if(comp[i] == -1)
			{
				comp[i]=nr;
				ComponentsR(nr, i, transposedGraph, comp);
			}
		}
	}
}


void DirectedGraph::PageRankA()
{
	std::random_device randDev;
	std::mt19937 gen(randDev());
	std::uniform_int_distribution<> chance(1, 100);	
	std::uniform_int_distribution<> start(0, m_numberOfVertices - 1);	

	int visitedNumber = 0;
	int nextStep;
	int currentVertice = start(gen);
	int randomNeighbourChoice;
	int visited = start(gen);
	double d = 15;
	double result[m_numberOfVertices];

	for (int i = 0; i < m_numberOfVertices; ++i)
	{
		result[i] = 0.0;
	}

	for(unsigned counter = 0; counter < 10000; counter++)
	{
		nextStep = chance(gen);
		if(nextStep <= (100 - d))
		{
			bool hasNeighbour = false;
			for (unsigned i = 0; i < m_numberOfVertices; ++i)
			{
				if(m_encodedGraphData.at(currentVertice).at(i) == 1)
				{
					hasNeighbour = true;
				}
			}	

			if(hasNeighbour)
			{
				hasNeighbour = false;
				do
				{
					randomNeighbourChoice = start(gen);
					if(m_encodedGraphData.at(currentVertice).at(randomNeighbourChoice) == 1)
					{
						result[randomNeighbourChoice] += 1;
						visitedNumber += 1;
						hasNeighbour = true;
						currentVertice = randomNeighbourChoice;
					}
				} while (!hasNeighbour);
			}

			else
			{
				randomNeighbourChoice = start(gen);
				result[randomNeighbourChoice] += 1;
				visitedNumber += 1;
				currentVertice = randomNeighbourChoice;
			}
		}

		else
		{
			randomNeighbourChoice = start(gen);
			result[randomNeighbourChoice] += 1;
			visitedNumber += 1;
			currentVertice = randomNeighbourChoice;
		}
	}

	std::cout<<"\nFrequency"<<std::endl;
	for (int i = 0; i < m_numberOfVertices; ++i)
	{
		if(i != m_numberOfVertices - 1)
			std::cout<<result[i]/visitedNumber<<", ";
		else
			std::cout<<result[i]/visitedNumber<<"\n";
	}	
}

void DirectedGraph::PageRankB()
{
	int n = m_numberOfVertices;
	int neighboursNumber;
	int di[n];
	double sum;
	double d = 0.15;
	double p[n];
	double p_temp[n];
	double P[n][n];

	for (int i = 0; i < n; ++i)
		{
			p[i] = 1/(double)n;
			//p_temp[i] = 1/(double)n;
		}

	for (int i = 0; i < n; ++i)
	{
		neighboursNumber = 0;
		for (int j = 0; j < n; ++j)
		{
			if(m_encodedGraphData.at(i).at(j) == 1)
				neighboursNumber += 1;
		}
		di[i] = neighboursNumber;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			P[i][j] = (1.0-d)*m_encodedGraphData.at(i).at(j)/(double)di[i] + d/(double)n;
		}
	}

	for(int t = 1; t < 1000; t++)
	{	
		for (int j = 0; j < n; ++j)
		{
			sum = 0;
			for (int k = 0; k < n; ++k)
			{
				sum += p[k] * P[k][j];
			}
			p_temp[j] = sum;
		}
		for (int i = 0; i < n; ++i)
			{
				p[i] = p_temp[i];
			}		
	}

	std::cout<<"\nFrequency\n";
	for (int i = 0; i < n; ++i)
	{			
			std::cout<<p[i]<< " ";		
	}
}

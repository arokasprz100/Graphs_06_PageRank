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
	m_hasWeights = 0;

	InitializeDistancesAndP();
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
	graph->InitializeDistancesAndP();

	return graph;
}

DirectedGraph* DirectedGraph::GenerateRandomStronglyConnectedGraph()
{
	std::random_device randDev;
	std::mt19937 gen(randDev());
	std::uniform_int_distribution<> ver(1, 10);
	std::uniform_int_distribution<> prob(1, 100);

	DirectedGraph * graph = DirectedGraph::GenerateRandomGraphBasedOnProbability(5, 1);

	while (graph->m_isStronglyConected == false)
	{
		delete graph;
		graph = DirectedGraph::GenerateRandomGraphBasedOnProbability(ver(gen), prob(gen));
		graph->Kosaraju();
	}
	return graph;
}

void DirectedGraph::InitializeDistancesAndP()
{
	for (unsigned i = 0; i < m_numberOfVertices; ++i) {
		p.emplace_back();
		distances.emplace_back();
		for (unsigned j = 0; j < m_numberOfVertices; ++j)
		{
			p[i].push_back(-1);
			distances[i].push_back(9999);
		}
	}
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

void DirectedGraph::PrintWeightMatrix(std::ostream& placeToPrint) const
{
	if (!m_hasWeights)
	{
		std::cout << "This graph has no weights, so they can not be printed" << std::endl;
		return;
	}

	placeToPrint << std::endl;
	placeToPrint << "Weight matrix" << std::endl << std::endl;

	for (unsigned i = 0; i<m_numberOfVertices; ++i)
	{
		for (unsigned j = 0; j< m_numberOfVertices; ++j)
			placeToPrint << m_weightMatrix.at(i).at(j) << " ";
		placeToPrint << std::endl;
	}
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



void DirectedGraph::GenerateWeights(int lowerBound, int upperBound)
{
	if (m_hasWeights)
	{
		std::cout<<"This graph already has weights. "<<std::endl;
		return;
	}

	std::random_device randDev;  
    std::mt19937 gen(randDev()); 
    std::uniform_int_distribution<> dis(lowerBound,upperBound);

    GraphRepresentation oldRepresentation = m_currentRepresentation;
    ChangeToAdjacencyMatrix();

    for (unsigned i =0; i<m_numberOfVertices; ++i)
    {
    	m_weightMatrix.emplace_back();
    	for (unsigned j =0; j<m_numberOfVertices; ++j)
    	{
    		if (m_encodedGraphData.at(i).at(j))
    			m_weightMatrix.at(i).push_back(dis(gen));
    		else
    			m_weightMatrix.at(i).push_back(0);
    	}
    }

    m_hasWeights = true;

    if (oldRepresentation == GraphRepresentation::AdjacencyList)
    	ChangeToAdjacencyList();
    else if (oldRepresentation == GraphRepresentation::IncidenceMatrix)
    	ChangeToIncidenceMatrix();
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
	
	// if(nr == 1)
		m_isStronglyConected = true;
	
	// std::cout<<"This graph has "<<nr<<" consistent components."<<std::endl;
	// for(int i=1; i<=nr; i++)
	// {
	// 	std::cout<<i<<" consistent component: ";
	// 	for(unsigned j = 0; j < comp.size(); ++j)
	// 	{
	// 		if(comp[j] == i)
	// 			std::cout<<j<<" ";
	// 	}
	// 	std::cout<<std::endl;
	// }

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

bool DirectedGraph::BellmanFord(int source, bool print)
{
	ChangeToAdjacencyMatrix();
	bool negativeCycle = false;

	Init(source);
	
	for (unsigned i = 1; i < m_numberOfVertices; ++i)
	{
		for (unsigned w = 0; w < m_numberOfVertices; ++w)
		{
			for (unsigned k = 0; k < m_numberOfVertices; ++k)
			{
				if(m_encodedGraphData[w][k]!=0)
					Relax(m_weightMatrix, w, k, source); 
			}
		}
	}

	for(unsigned w = 0; w < m_numberOfVertices; ++w)
	{
		for(unsigned k = 0; k < m_numberOfVertices; ++k)
			if((m_encodedGraphData[w][k]!=0) && (distances[source][k] > (distances[source][w]+m_weightMatrix[w][k])))
				negativeCycle = true;
	}

	if (negativeCycle)
	{
		std::cout << "Graph contains a negative cycle!" << std::endl;
		return false;
	}

	if (print == false)
		return true;

	std::cout<<"Vertex nr\tDistance"<<std::endl;
	for (unsigned i = 0; i < m_numberOfVertices; i++)
		std::cout<<i<<"\t\t"<<distances[source][i]<<std::endl;

	return true;
}



DirectedGraph* DirectedGraph::Add_S()
{
	DirectedGraph* extendedGraph = new DirectedGraph(*this);
	extendedGraph->ChangeToAdjacencyList();

	extendedGraph->m_encodedGraphData.emplace_back();
	extendedGraph->m_weightMatrix.emplace_back();
	extendedGraph->distances.emplace_back();
	extendedGraph->p.emplace_back();

	for (unsigned i = 0; i < m_numberOfVertices; ++i) 
	{
		extendedGraph->m_weightMatrix.at(m_numberOfVertices).push_back(0);
		extendedGraph->m_encodedGraphData.at(m_numberOfVertices).push_back(i);
		extendedGraph->m_weightMatrix.at(i).push_back(0);

		extendedGraph->distances[i].push_back(9999);
		extendedGraph->distances[m_numberOfVertices].push_back(9999);

		extendedGraph->p[i].push_back(-1);
		extendedGraph->p[m_numberOfVertices].push_back(-1);
	}

	
	extendedGraph->distances[m_numberOfVertices].push_back(9999);
	extendedGraph->p[m_numberOfVertices].push_back(-1);

	extendedGraph->m_weightMatrix.at(m_numberOfVertices).push_back(0);

	extendedGraph->m_numberOfVertices++;

	return extendedGraph;
}


void DirectedGraph::Init(int source)
{
	for (unsigned i = 0; i < m_numberOfVertices; ++i)
	{
		distances[source][i] = 9999;
		p[source][i] = -1;
	}
	distances[source][source] = 0;
}

void DirectedGraph::Relax(matrix2d& weights, int u, int v, int source)
{
	if (distances[source][v] > distances[source][u] + weights[u][v]) 
	{
		distances[source][v] = distances[source][u] + weights[u][v];
		p[source][v] = u;
	}
}

void DirectedGraph::Dijkstra(matrix2d& weights, int source)
{
	Init(source);
	std::set<int> S;
	while (S.size() != m_numberOfVertices)
	{
		int min = 9999;
		for (unsigned i = 0; i < m_numberOfVertices; ++i)
		{
			if (distances[source][i] < min && S.find(i) == S.end()) {
				S.insert(i);
				ChangeToAdjacencyList();
				for (unsigned j = 0; j < m_encodedGraphData.at(i).size(); ++j)
					Relax(weights, i, m_encodedGraphData.at(i).at(j), source);
			}
		}
	}
}

matrix2d DirectedGraph::Johnson()
{
	DirectedGraph* gPrim = Add_S();
	int s = gPrim->m_numberOfVertices - 1;
	for (unsigned i = 0; i < gPrim->m_numberOfVertices; ++i)
		gPrim->distances[s][i] = 9999;
	gPrim->distances[s][m_numberOfVertices - 1] = 0;

	matrix2d wPrim = gPrim->m_weightMatrix;
	matrix2d D;

	for (unsigned i = 0; i < m_numberOfVertices; ++i)
	{
		D.emplace_back();
		for (unsigned j = 0; j < m_numberOfVertices; ++j)
			D[i].push_back(0);
	}

	std::vector <int> h;
	for (unsigned i = 0; i < gPrim->m_numberOfVertices; ++i)
		h.push_back(9999);

	if (gPrim->BellmanFord(s, false) == false)
	{
		std::cout << "G contains negative cycle." << std::endl;
		throw std::invalid_argument("WRONG GRAPH. G contains negative cycle.");
	}
	else
	{
		for (unsigned i = 0; i < m_numberOfVertices; ++i)
			h[i] = gPrim->distances[s][i];

		gPrim->ChangeToAdjacencyMatrix();

		for (unsigned i = 0; i < gPrim->m_numberOfVertices; ++i)
		{
			for (unsigned j = 0; j < gPrim->m_numberOfVertices; ++j)
			{
				if (gPrim->m_encodedGraphData[i][j])
					wPrim[i][j] = gPrim->m_weightMatrix[i][j] + h[i] - h[j];
			}
			for (unsigned j = 0; j < m_numberOfVertices; ++j)
			{
				Dijkstra(wPrim, j);
				for (unsigned k = 0; k < m_numberOfVertices; ++k)
				{
					D[j][k] = distances[j][k] - h[j] + h[k];
				}
			}
		}
	}

	std::cout << std::endl;
	std::cout << "D matrix: " << std::endl;
	for (unsigned i = 0; i < m_numberOfVertices; ++i) {
		for (unsigned j = 0; j < m_numberOfVertices; ++j)
			std::cout << D[i][j] << " ";
		std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << "p matrix:" << std::endl;
	for (unsigned i = 0; i < m_numberOfVertices; ++i)
	{
		for (unsigned j = 0; j < m_numberOfVertices; ++j)
		{
			std::cout << p[i][j] << " ";
		}
		std::cout << std::endl;
	}

	return D;
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

}

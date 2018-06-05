#include "DirectedGraph.h"
#include "FileInput.h"
#include <iostream>
#include <cmath>
#include <random>

int main()
{

	std::cout<<std::endl<<std::endl;
	std::cout<<"In order to test PageRank algorithm, load graph from file or choose to randomly generate it: "<<std::endl;
	std::cout<<"1. Load graph from file. "<<std::endl;
	std::cout<<"2. Generate random graph. "<<std::endl;

	DirectedGraph* graph = nullptr;
	int userChoice;

	if (!(std::cin>>userChoice))
		{
			std::cin.clear();
			std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

	if(userChoice == 1)
	{
		if (graph)
			delete graph;
		graph = LoadGraphFromFile();
		graph->PrintGraph(std::cout);

		std::cout<<"\nExercise a)";
		graph->PageRankA();

		std::cout<<"\nExercise b)";
		graph->PageRankB();
	}

	if (userChoice == 2)
	{
		
		if (graph)
			delete graph;

		graph = DirectedGraph::GenerateRandomStronglyConnectedGraph();
		graph->PrintGraph(std::cout);

		std::cout<<"\nExercise a)";
		graph->PageRankA();

		std::cout<<"\nExercise b)";
		graph->PageRankB();
	}
	std::cout<<std::endl;
	if (graph)
		delete graph;

	return 0;
}

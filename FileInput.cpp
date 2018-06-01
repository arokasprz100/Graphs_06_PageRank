#include "FileInput.h"
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>

DirectedGraph* LoadGraphFromFile()
{
	std::ifstream file;
	for (;;){
		if (OpenFile(GetFileNameFromUser(), file))
			break;
	}

	std::string lineFromFile;
	std::string reprTypeFromFile;

	std::getline(file, lineFromFile);
	reprTypeFromFile = lineFromFile;
	reprTypeFromFile.erase( std::remove(reprTypeFromFile.begin(), reprTypeFromFile.end(), '\r'), reprTypeFromFile.end() );

	if (!CheckIfInputFileValid(reprTypeFromFile))
		throw std::runtime_error("Wrong input file format. Aborting.");

	std::getline(file, lineFromFile); 

	matrix2d dataFromFile;
	unsigned currentVertex = 0;

	while (std::getline(file, lineFromFile))
	{
		file.clear();
		dataFromFile.emplace_back();

		int adjacentVertex;
		std::stringstream streamFromFile(lineFromFile);
		while (streamFromFile>>adjacentVertex)
			dataFromFile.at(currentVertex).push_back(adjacentVertex);

		currentVertex++;

	}

	file.close();

	return new DirectedGraph(dataFromFile, reprTypeFromFile);
}

std::string GetFileNameFromUser()
{
	std::string fileName;
	std::cout<<"Please enter proper file name: ";
	std::getline(std::cin, fileName);
	return fileName;
}

bool OpenFile(const std::string& fileName, std::ifstream& file)
{
	file.open(fileName.c_str());
	if (file.is_open())
		return true;
	else return false;
}

bool CheckIfInputFileValid(std::string firstLine)
{
	if (firstLine == "adjacency_matrix")
		return true;
	else if (firstLine == "adjacency_list")
		return true;
	else if (firstLine == "incidence_matrix")
		return true;
	else return false;
}
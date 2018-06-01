#ifndef FileInput_h
#define FileInput_h

#include "DirectedGraph.h"
#include <iostream>
#include <fstream>

DirectedGraph* LoadGraphFromFile();

std::string GetFileNameFromUser();

bool OpenFile(const std::string& fileName, std::ifstream& file);

bool CheckIfInputFileValid(std::string firstLine);

#endif
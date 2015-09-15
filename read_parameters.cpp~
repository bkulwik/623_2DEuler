
//This is used to read all parameters from the input file
#include "file_header.h"
#include <iostream>
#include <fstream>

std::vector<double> read_parameters(std::string filename) {
	std::vector<double> inputs;
	std::string line;

// Read input file with defined solver parameters
// First - CFL, second - tmax
	std::ifstream parameter_file(filename);
	if (parameter_file.is_open()) {
		while (getline (parameter_file,line)) {
			if (line[0] != '%') {
				inputs.push_back(std::stod(line));
			}
		}
	parameter_file.close();
	}
	else {
		std::cout << "Cannot Open File" << '\n';
	}
return inputs;
}


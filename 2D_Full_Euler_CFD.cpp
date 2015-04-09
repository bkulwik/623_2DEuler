#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "file_header.h"

int main() {
// Define Variables
	std::string parameter_filename = "solver_parms.txt";
	std::string input_filename = "test_input.bkcfd";
	std::vector<double> parameters;
	double CFL;
	double tmax;

	std::vector<std::vector<int>> cell_matrix;
	std::vector<std::vector<double>> node_matrix;
	std::vector<std::vector<double>> U_matrix;
	
	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	CFL = parameters[0];	
	tmax = parameters[1];
	
/*
	std::cout << "Parameters are: " << '\n';
	for(auto input:parameters) {
		std::cout << input << '\n';
	}
	std::cout << '\n';
*/

// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v e] defined for each cell


	read_grid(input_filename, node_matrix, cell_matrix, U_matrix);

	for(auto input:node_matrix) {
		std::cout << input[0] << ' ';
		std::cout << input[1] << ' ';
		std::cout << input[2] << '\n';
	}



// Start calculation in for loop going to final time

//for (int t = 0, t < tmax; t++) {
 
// Go through each edge, calculate flux, save

// for (int edge = 1, edge < numedges, edge++) {

// Solve riemann problem on edge for euler flux


//}

//Use all these fluxes to define updated state on each cell 

//Option to save at each timestep


//}
return 0;
}





 

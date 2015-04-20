// This is used to read in all solver parameters in the solver_parms.txt file
// It may need to be updated with other parameters as we go, and as other parameters need to be taken into account in this code.


#ifndef FILE_HEADER_H
#define FILE_HEADER_H
#include <vector>
#include <string>

struct ODstate { // ODstate for 1-dimensional state
	double rho;
	double rhou;
	double E;
	double pressure;
};

struct TDstate { // TDstate for 2-dimensional state
	double rho;
	double rhou;
	double rhov;
	double E;
	double pressure;
};

struct cell {
	int cellnumber;
	std::vector<double> cornerlocs_x;
	std::vector<double> cornerlocs_y;
	std::vector<int> adjacent_cells;
	TDstate state;
};


//This is used to read the parameters from the parameter file
std::vector<double> read_parameters(std::string filename);

//This is used to read the input node locations and cell connectivity from the given input_filename
void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma);

// This is the 1D exact riemann solver
ODstate Exact_Riemann_Solver(ODstate left, ODstate right, double thresh, double gamma);

// Given a cell number, this computes the cells adjacent to it
//cell_neighbors compute_neighbors(int cell, std::vector<std::vector<int>> cell_matrix, std::vector<std::vector<double>> node_matrix);

#endif
//==============================================================

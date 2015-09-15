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
	int edge_type;
	TDstate state;
};


//This is used to read the parameters from the parameter file
std::vector<double> read_parameters(std::string filename);

//This is used to read the input node locations and cell connectivity from the given input_filename
void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma);

// This is the 1D exact riemann solver
ODstate Exact_Riemann_Solver(ODstate left, ODstate right, double left_parallelvel, double right_parallelvel, double thresh, double gamma, bool debug);

//This is used to write solver outputs to a .txt file to be read into and plotted by MATLAB
void write_to_file(const std::vector<cell>& grid, const std::vector<TDstate>& U, std::string& filename);

#endif
//==============================================================

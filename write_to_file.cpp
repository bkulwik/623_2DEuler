#include "file_header.h"
#include "helper_functions.h"
#include <iostream>
#include <fstream>


void output_cornerlocs(const std::vector<double>& cornerlocs, std::ofstream& ofile);

void output_state(const TDstate& U, std::ofstream& ofile);

void write_to_file(const std::vector<cell>& grid, const std::vector<TDstate>& U, std::string& filename) {
	std::ofstream ofile;
	ofile.open(filename);
	
	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) {
		output_cornerlocs(grid[cellnum].cornerlocs_x, ofile);
		output_cornerlocs(grid[cellnum].cornerlocs_y, ofile);
		output_state(U[cellnum], ofile);
		ofile << '\n';
	}
	ofile.close();
}

void output_cornerlocs(const std::vector<double>& cornerlocs, std::ofstream& ofile) {
	for(unsigned int cornernum = 0; cornernum < 4; ++cornernum) {
		ofile << cornerlocs[cornernum] << ", ";
	}
}

void output_state(const TDstate& U, std::ofstream& ofile) {
	ofile << U.rho << ", " << U.rhou << ", " << U.rhov << ", " << U.E;
}


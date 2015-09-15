#include <iostream>
#include <fstream>
#include <sstream>
#include "file_header.h"
#include "helper_functions.h"

// Pass by reference: std::vector<int>& cell_number

// This reads in the grid from the input file specified in input_filename
//
// Input file should be in this format:
//
// "cellnumber","cornerlocs_x","cornerlocs_y","adjacent_cells","initial_conditions",
// 0,0,0.125,0.125,0,0,0,1.3333,1.3333,0,0,2,9,1,0,0,250000,
/*
const std::string input_file = "firstline\n\
0,0,1,1,0,0,0,1,1,-1,-1,1,5,1,1,0,0,200000,\n\
1,1,2,2,1,0,0,1,1,0,-1,2,6,1,1,0,0,200000,\n\
2,2,3,3,2,0,0,1,1,1,-1,3,7,1,1,0,0,200000,\n\
3,3,4,4,3,0,0,1,1,2,-1,4,8,1,1,0,0,200000,\n\
4,4,5,5,4,0,0,1,1,3,-1,-1,9,1,1,0,0,200000,\n\
5,0,1,1,0,1,1,2,2,-1,0,6,10,1,1,0,0,200000,\n\
6,1,2,2,1,1,1,2,2,5,1,7,11,0,1,0,0,200000,\n\
7,2,3,3,2,1,1,2,2,6,2,8,12,0,1,0,0,200000,\n\
8,3,4,4,3,1,1,2,2,7,3,9,13,0,1,0,0,200000,\n\
9,4,5,5,4,1,1,2,2,8,4,-1,14,1,1,0,0,200000,\n\
10,0,1,1,0,2,2,3,3,-1,5,11,15,1,1,0,0,200000,\n\
11,1,2,2,1,2,2,3,3,10,6,12,16,0,1,0,0,200000,\n\
12,2,3,3,2,2,2,3,3,11,7,13,17,0,1,0,0,200000,\n\
13,3,4,4,3,2,2,3,3,12,8,14,18,0,1,0,0,200000,\n\
14,4,5,5,4,2,2,3,3,13,9,-1,19,1,1,0,0,200000,\n\
15,0,1,1,0,3,3,4,4,-1,10,16,20,1,1,0,0,200000,\n\
16,1,2,2,1,3,3,4,4,15,11,17,21,0,1,0,0,200000,\n\
17,2,3,3,2,3,3,4,4,16,12,18,22,0,1,0,0,200000,\n\
18,3,4,4,3,3,3,4,4,17,13,19,23,0,1,0,0,200000,\n\
19,4,5,5,4,3,3,4,4,18,14,-1,24,1,1,0,0,200000,\n\
20,0,1,1,0,4,4,5,5,-1,15,21,-1,1,1,0,0,200000,\n\
21,1,2,2,1,4,4,5,5,20,16,22,-1,1,1,0,0,200000,\n\
22,2,3,3,2,4,4,5,5,21,17,23,-1,1,1,0,0,200000,\n\
23,3,4,4,3,4,4,5,5,22,18,24,-1,1,1,0,0,200000,\n\
24,4,5,5,4,4,4,5,5,23,19,-1,-1,1,1,0,0,200000,";
*/

void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma)  {
	std::string line;
	std::cout << "reading grid... \n";
	int sep_loc_1;
	int sep_loc_old, sep_loc_new;
	int linenumber = 1;
	int cell_num = 0;

	std::ifstream input_file(input_filename);
	if (input_file.is_open()) {
		while (getline(input_file,line)) {
			if (linenumber == 1) {
				linenumber++;
				continue;
			}
			else if (!line.empty()) {
				grid.push_back(cell());
				// Read cell number
				sep_loc_1 = line.find(',');
				std::string cell_number_str = line.substr(0, sep_loc_1);
				grid[cell_num].cellnumber = (std::stoi(cell_number_str));

				// Read corner node x locations
				sep_loc_old = sep_loc_1;
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string node_x_temp = line.substr(sep_loc_old+1, sep_loc_new);		
					grid[cell_num].cornerlocs_x.push_back(std::stod(node_x_temp));
					sep_loc_old = sep_loc_new;
				}

				// Read corner node y locations				
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string node_y_temp = line.substr(sep_loc_old+1, sep_loc_new);		
					grid[cell_num].cornerlocs_y.push_back(std::stod(node_y_temp));
					sep_loc_old = sep_loc_new;
				}
				
				// Read adjacent cells				
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string adjacentcell_temp = line.substr(sep_loc_old+1, sep_loc_new);		
					grid[cell_num].adjacent_cells.push_back(std::stod(adjacentcell_temp));
					sep_loc_old = sep_loc_new;
				}

				// Read edge_cell boolean
				sep_loc_new = line.find(',',sep_loc_old+1);
				std::string adjacentcell_temp = line.substr(sep_loc_old+1, sep_loc_new);		
				grid[cell_num].edge_type = (std::stoi(adjacentcell_temp));
				sep_loc_old = sep_loc_new;

				// Read initial conditions				
				for(int iX = 1; iX <= 4; ++iX) {
					sep_loc_new = line.find(',',sep_loc_old+1);
					std::string initialcondition_temp = line.substr(sep_loc_old+1, sep_loc_new);		
					if (iX == 1) {
						grid[cell_num].state.rho = (std::stod(initialcondition_temp));
					} else if (iX == 2) {
						grid[cell_num].state.rhou = (std::stod(initialcondition_temp));
					} else if (iX == 3) {
						grid[cell_num].state.rhov = (std::stod(initialcondition_temp));
					} else if (iX == 4) {
						grid[cell_num].state.E = (std::stod(initialcondition_temp));
					}
					sep_loc_old = sep_loc_new;
				}
				grid[cell_num].state.pressure = compute_pressure_2D(grid[cell_num].state, gamma);

			}
			else {
				break;
			}
			linenumber++;
			cell_num++;
		}
	} else {
		std::cout << "Input file not open... file does not exist." << '\n';
	}
}

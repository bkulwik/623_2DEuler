#include <iostream>
#include <fstream>
#include "file_header.h"
#include "helper_functions.h"

// Pass by reference: std::vector<int>& cell_number

// This reads in the grid from the input file specified in input_filename
//
// Input file should be in this format:
//
// "cellnumber","cornerlocs_x","cornerlocs_y","adjacent_cells","initial_conditions",
// 0,0,0.125,0.125,0,0,0,1.3333,1.3333,0,0,2,9,1,0,0,250000,


void read_grid(std::string input_filename, std::vector<cell> &grid, double gamma)  {
	std::string line;
	
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
	}
}
/*
				sep_loc_3 = line.find(' ',sep_loc_2+1);
				std::string node_y_temp = line.substr(sep_loc_2+1, sep_loc_3);
				node_y = (std::stod(node_y_temp));
				node_temp.at(2) = node_y;

				node_matrix.push_back(node_temp);
			}
			else if ((line[0] != '%') && (entered == 1)) {
			// Read cell numbers and counterclockwise node numbers, starting from bottom left corner of cell	
				if (line.empty()) {
					// Exit here and start looking for initial conditions
					entered = 2;
					break;
				}

				sep_loc_1 = line.find(' ');
				std::string cell_number_temp = line.substr(0, sep_loc_1);		
				cell_number = (std::stod(cell_number_temp));
				cell_temp.at(0) = cell_number;

				sep_loc_2 = line.find(' ',sep_loc_1+1);
				std::string cell_BL_temp = line.substr(sep_loc_1+1, sep_loc_2);
				cell_BL = (std::stod(cell_BL_temp));
				cell_temp.at(1) = cell_BL;

				sep_loc_3 = line.find(' ',sep_loc_2+1);
				std::string cell_BR_temp = line.substr(sep_loc_2+1, sep_loc_3);
				cell_BR = (std::stod(cell_BR_temp));
				cell_temp.at(2) = cell_BR;

				sep_loc_4 = line.find(' ',sep_loc_3+1);
				std::string cell_TR_temp = line.substr(sep_loc_3+1, sep_loc_4);
				cell_TR = (std::stod(cell_TR_temp));
				cell_temp.at(3) = cell_TR;

				sep_loc_5 = line.find(' ',sep_loc_4+1);
				std::string cell_TL_temp = line.substr(sep_loc_4+1, sep_loc_5);
				cell_TL = (std::stod(cell_TL_temp));
				cell_temp.at(4) = cell_TL;

				cell_matrix.push_back(cell_temp);
			}
			else if ((line[0] != '%') && (entered == 2)) {
				if (line.empty()) {
					//We are done with this file
					break;
				}

				sep_loc_1 = line.find(' ');
				std::string cell_number_temp = line.substr(0, sep_loc_1);		
				cell_number = (std::stod(cell_number_temp));
				U_temp.at(0) = cell_number;

				sep_loc_2 = line.find(' ',sep_loc_1+1);
				std::string rho_temp = line.substr(sep_loc_1+1, sep_loc_2);
				rho = (std::stod(rho_temp));
				U_temp.at(1) = rho;

				sep_loc_3 = line.find(' ',sep_loc_2+1);
				std::string rhou_temp = line.substr(sep_loc_2+1, sep_loc_3);
				rhou = (std::stod(rhou_temp));
				U_temp.at(2) = rhou;

				sep_loc_4 = line.find(' ',sep_loc_3+1);
				std::string rhov_temp = line.substr(sep_loc_3+1, sep_loc_4);
				rhov = (std::stod(rhov_temp));
				U_temp.at(3) = rhov;

				sep_loc_5 = line.find(' ',sep_loc_4+1);
				std::string e_temp = line.substr(sep_loc_4+1, sep_loc_5);
				e = (std::stod(e_temp));
				U_temp.at(4) = e;

				U_matrix.push_back(U_temp);

			}
		}
	}
*/
/*
std::cout << node_matrix[0][0] << ' ' << node_matrix[0][1] << ' ' << node_matrix[0][2] << '\n';
std::cout << node_matrix[1][0] << ' ' << node_matrix[1][1] << ' ' << node_matrix[1][2] << '\n';
std::cout << node_matrix[2][0] << ' ' << node_matrix[2][1] << ' ' << node_matrix[2][2] << '\n';
std::cout << node_matrix[3][0] << ' ' << node_matrix[3][1] << ' ' << node_matrix[3][2] << '\n';
std::cout << node_matrix[4][0] << ' ' << node_matrix[4][1] << ' ' << node_matrix[4][2] << '\n';
std::cout << node_matrix[5][0] << ' ' << node_matrix[5][1] << ' ' << node_matrix[5][2] << '\n';
std::cout << node_matrix[6][0] << ' ' << node_matrix[6][1] << ' ' << node_matrix[6][2] << '\n';
std::cout << node_matrix[7][0] << ' ' << node_matrix[7][1] << ' ' << node_matrix[7][2] << '\n';
std::cout << node_matrix[8][0] << ' ' << node_matrix[8][1] << ' ' << node_matrix[8][2] << '\n';

std::cout << '\n';
std::cout << cell_matrix[0][0] << ' ' << cell_matrix[0][1] << ' ' << cell_matrix[0][2] << ' ' << cell_matrix[0][3] << ' ' << cell_matrix[0][4] << '\n';
std::cout << cell_matrix[1][0] << ' ' << cell_matrix[1][1] << ' ' << cell_matrix[1][2] << ' ' << cell_matrix[1][3] << ' ' << cell_matrix[1][4] << '\n';
std::cout << cell_matrix[2][0] << ' ' << cell_matrix[2][1] << ' ' << cell_matrix[2][2] << ' ' << cell_matrix[2][3] << ' ' << cell_matrix[2][4] << '\n';
std::cout << cell_matrix[3][0] << ' ' << cell_matrix[3][1] << ' ' << cell_matrix[3][2] << ' ' << cell_matrix[3][3] << ' ' << cell_matrix[3][4] << '\n';
*/



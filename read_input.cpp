#include "file_header.h"
#include <iostream>
#include <fstream>


// Pass by reference: std::vector<int>& cell_number

void read_grid(std::string input_filename, std::vector<std::vector<double>> &node_matrix, std::vector<std::vector<int>> &cell_matrix std::vector<std::vector<double>> &U_matrix)  {
	
	std::string line;
	int entered = 0;
	
	int node_number;
	double node_x;
	double node_y;
	
	std::vector<double> node_temp(3);
	std::vector<int> cell_temp(5);

	int sep_loc_1;
	int sep_loc_2;
	int sep_loc_3;
	int sep_loc_4;
	int sep_loc_5;

	int cell_number;
	int cell_BL;
	int cell_BR;
	int cell_TR;
	int cell_TL;

	double rho
	double rhou
	double rhov
	double e
	
	std::ifstream input_file(input_filename);

	if (input_file.is_open()) {
		while (getline(input_file,line)) {
			if ((line[0] != '%') && (entered == 0)) {
				if (line.empty()) {
				// exit here and start looking for the cell definitions
					entered = 1;
					continue;
				}		
				// Read node numbers
				sep_loc_1 = line.find(' ');
				std::string node_number_temp = line.substr(0, sep_loc_1);		
				node_number = (std::stod(node_number_temp));
				node_temp.at(0) = node_number;


				// Read node x-y locations

				sep_loc_2 = line.find(' ',sep_loc_1+1);
				std::string node_x_temp = line.substr(sep_loc_1+1, sep_loc_2-2);		
				node_x = (std::stod(node_x_temp));
				node_temp.at(1) = node_x;

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

}

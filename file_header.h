// This is used to read in all solver parameters in the solver_parms.txt file
// It may need to be updated with other parameters as we go, and as other parameters need to be taken into account in this code.


#ifndef FILE_HEADER_H
#define FILE_HEADER_H
#include <vector>
#include <string>

//This is used to read the parameters from the parameter file
std::vector<double> read_parameters(std::string filename);

//This is used to read the input node locations and cell connectivity from the given input_filename
void read_grid(std::string input_filename, std::vector<std::vector<double>> &node_matrix, std::vector<std::vector<int>> &cell_matrix, std::vector<std::vector<double>> &U_matrix);

#endif
//==============================================================

#include "file_header.h"

// Computes 2D pressure from a TDstate
double compute_pressure_2D(TDstate state, double gamma);

// Computes 1D pressure from a ODstate
double compute_pressure_1D(ODstate state, double gamma);

// Computes the max value in a vector of doubles
double vectormax(std::vector<double> locations);

// Computes the min value in a vector of doubles
double vectormin(std::vector<double> locations);

// Computes the min delta_(x or y) value in grid
//double gridmin(std::vector<cell> grid, char direction);

// Computes and returns the cell numbers of all non-edge cells in the grid
std::vector<int> find_interior_cells(const std::vector<cell>& grid);

// Updates the max wavespeed from the given cell's left, right, top and bottom flux.
void max_wavespeed_calculator(double &max_wavespeed, std::vector<TDstate> &F, std::vector<TDstate> &G, TDstate state);

// Updates the max wavespeed from the riemann wavespeeds
void max_wavespeed_calculator_riemann(double& max_wavespeed, TDstate state, double gamma);

// Outputs the global state U to the terminal
void outputU(std::vector<TDstate>& U);

// Outputs a ODstate to the terminal
void outputODstate(ODstate& state);

// Outputs a TDstate to the terminal
void outputTDstate(TDstate& state);



//==============================================================

#include "file_header.h"

// Computes 2D pressure 
double compute_pressure_2D(TDstate state, double gamma);

// Computes 1D pressure
double compute_pressure_1D(ODstate state, double gamma);

// Computes the max value in a vector of doubles
double vectormax(std::vector<double> locations);

// Computes the min value in a vector of doubles
double vectormin(std::vector<double> locations);

// Computes the min delta_(x or y) value in grid
double gridmin(std::vector<cell> grid, char direction);

// Updates the max wavespeed from the given cell's left, right, top and bottom flux.
void max_wavespeed_calculator(double &max_wavespeed, std::vector<TDstate> &F, std::vector<TDstate> &G, TDstate state);

// Updates the max wavespeed from the riemann wavespeeds
void max_wavespeed_calculator_riemann(double& max_wavespeed, TDstate state, double gamma);

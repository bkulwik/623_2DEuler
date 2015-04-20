#include <iostream> 
#include <vector>
#include <cmath>
#include <cassert>
#include "helper_functions.h"

double compute_pressure_2D(TDstate state, double gamma) { 
// This computes the pressure for a 2D state
	double vel_mag = sqrt((pow(state.rhou,2) + pow(state.rhov,2))/(pow(state.rho,2)));
	double pressure = (state.E + state.rho*pow(vel_mag,2)/2)*(gamma-1);

	return(pressure);
}

double compute_pressure_1D(ODstate state, double gamma) {
	double pressure = (state.E + pow(state.rhou,2)/(2*state.rho))*(gamma-1);
	return(pressure);
}

double vectormax(std::vector<double> locations) {
	double max = -pow(10,50);
	for (unsigned int it = 0; it < locations.size(); it++) {
		if (it == 0) {
			max = locations[it];
		} 
		else if (locations[it] > max) {
			max = locations[it];
		}
	}
	return(max);
}

double vectormin(std::vector<double> locations) {
	double min_value = pow(10,50);
	for (unsigned int it = 0; it < locations.size(); it++) {
		if (it == 0) {
			min_value = locations[it];
		} 
		else if (locations[it] < min_value) {
			min_value = locations[it];
		}
	}
	return(min_value);
}

double gridmin(std::vector<cell> grid, char direction) {
	double delta, min_delta;
	if (direction == 'x') {
		min_delta = std::abs(grid[0].cornerlocs_x[0] - grid[0].cornerlocs_x[1]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			delta = std::abs(grid[it].cornerlocs_x[0] - grid[it].cornerlocs_x[1]);
			if (delta < min_delta) {
				min_delta = delta;
			}
			delta = std::abs(grid[it].cornerlocs_x[3] - grid[it].cornerlocs_x[2]);
			if (delta < min_delta) {
				min_delta = delta;
			}
		}
	} else { // direction == "y"
		min_delta = std::abs(grid[0].cornerlocs_y[3] - grid[0].cornerlocs_y[0]);
		for (unsigned int it = 0; it < grid.size(); ++it) {
			delta = std::abs(grid[it].cornerlocs_y[3] - grid[it].cornerlocs_y[0]);
			if (delta < min_delta) {
				min_delta = delta;
			}
			delta = std::abs(grid[it].cornerlocs_y[2] - grid[it].cornerlocs_y[1]);
			if (delta < min_delta) {
				min_delta = delta;
			}
		}
	}
	return(min_delta);
}

void max_wavespeed_calculator(double &max_wavespeed, std::vector<TDstate> &F, std::vector<TDstate> &G, TDstate state) {
	assert(F.size() == 2);
	assert(G.size() == 2);

	double wavespeed_temp;
	for (int it = 0; it < 2; ++it) {
		wavespeed_temp = (F[it].rho/state.rho);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 1: " << max_wavespeed << '\n';
		wavespeed_temp = (G[it].rho/state.rho);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 2: " << max_wavespeed << '\n';
		wavespeed_temp = (F[it].rhou/state.rhou);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
		std::cout << "Flux: " << F[it].rhou << ", State: " << state.rhou << '\n';
	std::cout << "temp max wavespeed 3: " << max_wavespeed << '\n';
		wavespeed_temp = (G[it].rhou/state.rhou);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 4: " << max_wavespeed << '\n';
		wavespeed_temp = (F[it].rhov/state.rhov);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 5: " << max_wavespeed << '\n';
		wavespeed_temp = (G[it].rhov/state.rhov);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}	
	std::cout << "temp max wavespeed 6: " << max_wavespeed << '\n';
		wavespeed_temp = (F[it].E/state.E);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 7: " << max_wavespeed << '\n';
		wavespeed_temp = (G[it].E/state.E);
		if (wavespeed_temp > max_wavespeed) {
			max_wavespeed = wavespeed_temp;
		}
	std::cout << "temp max wavespeed 8: " << max_wavespeed << '\n';
	}
}

void max_wavespeed_calculator_riemann(double& max_wavespeed, TDstate state, double gamma) {
	double a = sqrt(gamma*compute_pressure_2D(state, gamma)/state.rho);
	double U = state.rhou/state.rho;
	double V = state.rhov/state.rho;	
	std::vector<double> local_wavespeeds;
	local_wavespeeds.push_back(U+a);
	local_wavespeeds.push_back(U);
	local_wavespeeds.push_back(U-a);
	local_wavespeeds.push_back(V+a);
	local_wavespeeds.push_back(V);
	local_wavespeeds.push_back(V-a);

	double max_local_wavespeed = vectormax(local_wavespeeds);
	if (max_local_wavespeed > max_wavespeed) {
		max_wavespeed = max_local_wavespeed;
	}
}
	

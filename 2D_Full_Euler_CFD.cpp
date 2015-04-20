#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include "file_header.h"
#include "helper_functions.h"

void FluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellnum, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double min_delta, double CFL, double& max_wavespeed);

// Define helper functions!
bool isedge(const std::vector<cell>& grid, int cellnum);

// void compute_flux(std::vector<double> &fluxph, std::vector<double> &fluxmh, std::vector<double> slope, TDstate state, double delta, std::string direction, double gamma);

// This is not math-y, just a formality of grabbing the left/right (state_minus) or bottom/top (state_plus) states from the grid input.
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<cell>& grid, int cellnum, int direction);

// Computes the "geometric" slope vector for all conserved quantities between the bottom/left and top/right cells
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, std::vector<cell> grid, int cellnum, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y);

// This computes the double-value "limiter value" based on the bottom/left and top/right states for a given cell. The limiter value is then used to calculate the flux
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number);

// Computes the harmonic limiter given r, the ratio of minus to plus slopes
void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus);

void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, int direction, double gamma);

// Computes the value of flux in the center cell, not along the edges but in the cell
void compute_center_flux(TDstate &center_flux, TDstate center_state, int direction, double gamma);

// Reconstructs the solution from t to t+dt/2
void compute_halfway_state(unsigned int cellnum, TDstate& Uph, const TDstate& current_U, const std::vector<double>& F, const std::vector<double>& G, double dt, double dx, double dy);

int main() {

// Define Variables
	std::string parameter_filename = "solver_parms.txt";
	std::string input_filename = "5x5square2.bkcfd";
	std::vector<double> parameters;
	std::vector<double> Fiph, Fimh, Giph, Gimh;
	double CFL, tmax, dt, delta_x, delta_y, min_delta_x, min_delta_y, min_delta;
	TDstate state_minus, state_plus, slope_minus, slope_plus, limiter_value;
	std::vector<TDstate> F (2);
	std::vector<TDstate> G (2);

//	double thresh = 0.00001; // 1E-5 tolerance for riemann solver - velocity
	double gamma = 1.4;

// 1 = harmonic mean for slopes, 2 = ...
	int limiter_number = 1; 

	std::vector<cell> grid;

	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	CFL = parameters.at(0);	
	tmax = parameters.at(1);
	

// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v e] defined for each cell
	read_grid(input_filename, grid, gamma);
	
	min_delta_x = gridmin(grid, 'x');
	min_delta_y = gridmin(grid, 'y');
	min_delta = std::min(min_delta_x,min_delta_y);
	

/*		std::cout << "printing now..." << '\n';
	for(auto input: grid) {
		std::cout << input.cellnumber << '\n';
		std::cout << input.cornerlocs_x[0] << ' ';
		std::cout << input.cornerlocs_x[1] << ' ';
		std::cout << input.cornerlocs_x[2] << ' ';
		std::cout << input.cornerlocs_x[3] << '\n';
		std::cout << input.cornerlocs_y[0] << ' ';
		std::cout << input.cornerlocs_y[1] << ' ';
		std::cout << input.cornerlocs_y[2] << ' ';
		std::cout << input.cornerlocs_y[3] << '\n';
		std::cout << input.adjacent_cells[0] << ' ';
		std::cout << input.adjacent_cells[1] << ' ';
		std::cout << input.adjacent_cells[2] << ' ';
		std::cout << input.adjacent_cells[3] << '\n';
		std::cout << input.state.rho << ' ';
		std::cout << input.state.rhou << ' ';
		std::cout << input.state.rhov << ' ';
		std::cout << input.state.E << '\n' << '\n';
	} 
*/

// Start calculation in for loop going to final time
double t = 0;
//	for (double t = 0, t <= tmax; t++) 
	double max_wavespeed = 0; 
	
	if (t == 0) { // Compute initial dt for first timestep
		for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
			if (!isedge(grid,cellnum)) {
            FluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellnum, limiter_value, F, G, min_delta, CFL, max_wavespeed);
			}
		}
		dt = CFL*min_delta/max_wavespeed;
	}

// Go through each cell, calculate fluxes, update to new piecewise linear state
	std::vector<TDstate> Uph (grid.size());	

	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
		if (!isedge(grid,cellnum)) {
            FluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellnum, limiter_value, F, G, min_delta, CFL, max_wavespeed);
		// Now, use the fluxes calculated above to assign a new state, Uph
		compute_halfway_state(cellnum, Uph, grid[cellnum].state, F, G, dt, delta_x, delta_y);
									
		} // end of if (!isedge(grid,cell)) 	


// -----------------------------
// Need to check values of flux to see if they are giving accurate results!

//
// Go through each edge, calculate flux, save
//		for (int edge = 1, edge <= numedges, ++edge)
//
// Solve riemann problem on edge for euler flux
		ODstate left;
		ODstate right;
		left.rho = 1.0; //rho
		left.rhou = 300; //rho*u or rho*v
		left.E = 100000/0.4 - pow(left.rhou,2)/(left.rho*2); //E
		right.rho = 1;
		right.rhou = 0;
		right.E = 100000/0.4 - pow(right.rhou,2)/(right.rho*2);				

/*		ODstate test_state = Exact_Riemann_Solver(left, right, thresh, gamma);

		std::cout << '\n' << "State on x=0 is:" << '\n';
		std::cout << test_state.rho << '\n';
		std::cout << test_state.rhou << '\n';
		std::cout << test_state.E << '\n';
		std::cout << test_state.pressure << '\n';
*/


//Use all these fluxes to define updated state on each cell 

//Option to save at each timestep

	} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell)
return 0;
}


// This is the land of the helper functions! ---------------------------------------------------

void FluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellnum, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double min_delta, double CFL, double& max_wavespeed) {

	delta_x = (vectormax(grid[cellnum].cornerlocs_x) - vectormin(grid[cellnum].cornerlocs_x));
	delta_y = (vectormax(grid[cellnum].cornerlocs_y) - vectormin(grid[cellnum].cornerlocs_y));

    // Calculate the limited flux of each conserved quantity for initial t -> t+1/2 update
	for (int direction = 0; direction < 2; ++direction) { // For the left/right direction and up/down direction
		// Assign the bottom/left and the top/right states for each cell
		assign_edge_state(state_minus, state_plus, grid, cellnum, direction);	
						
		// Calculate the flux of each conserved variable, in each direction, on each face, of that cell.
		// Flux needs to be calculated by first finding the limiter (direction-dependent), then using that limiter to calculate the directional flux, then using the two directional fluxes to update from u_t to u_t+1/2	
			
		compute_slope(slope_minus, slope_plus, grid, cellnum, direction, state_minus, state_plus, delta_x, delta_y);
		compute_limiter(limiter_value, slope_minus, slope_plus, limiter_number);

		if (direction == 0) { // Compute F, x-fluxes
			compute_limited_flux(F, limiter_value, grid[cellnum].state, slope_plus, direction, gamma);
//			std::cout << "rho limiter x: " << limiter_value.rho << '\n';
		} else { // Compute G, y-fluxes
//			std::cout << "rho limiter y: " << limiter_value.rho << '\n';
			compute_limited_flux(G, limiter_value, grid[cellnum].state, slope_plus, direction, gamma);
		}
	} // end of for (int direction = 0; direction < 2; ++direction)
//	std::cout << "F flux: " << F[0].rho << " " << F[1].rho << ", G flux: " << G[0].rho << " " << G[1].rho << '\n';	
	max_wavespeed_calculator_riemann(max_wavespeed, grid[cellnum].state, gamma);
}


bool isedge(const std::vector<cell>& grid, int cellnum) {
	int special_value = -1;
	if ((grid[cellnum].adjacent_cells[0] == special_value) || (grid[cellnum].adjacent_cells[1] == special_value) || (grid[cellnum].adjacent_cells[2] == special_value)  || (grid[cellnum].adjacent_cells[3] == special_value)) {
		return(1);
	}
	else {
		return(0);
	}
}


void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<cell>& grid, int cellnum, int direction) {
	int adjacent_bl, adjacent_tr;

	adjacent_bl = (grid[cellnum].adjacent_cells[direction]);
	adjacent_tr = (grid[cellnum].adjacent_cells[direction+2]);

	state_minus = (grid[adjacent_bl].state);
	state_plus = (grid[adjacent_tr].state);
}


void compute_slope(TDstate &slope_minus, TDstate &slope_plus, std::vector<cell> grid, int cellnum, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y) {
	if (direction == 0) { //left-right
		slope_minus.rho = ((grid[cellnum].state.rho - state_minus.rho)/delta_x); // was 2*delta_x, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((grid[cellnum].state.rhou - state_minus.rhou)/delta_x);
		slope_minus.rhov = ((grid[cellnum].state.rhov - state_minus.rhov)/delta_x);
		slope_minus.E = ((grid[cellnum].state.E - state_minus.E)/delta_x);
		slope_plus.rho = ((state_plus.rho - grid[cellnum].state.rho)/delta_x);
		slope_plus.rhou = ((state_plus.rhou - grid[cellnum].state.rhou)/delta_x);
		slope_plus.rhov = ((state_plus.rhov - grid[cellnum].state.rhov)/delta_x);
		slope_plus.E = ((state_plus.E - grid[cellnum].state.E)/delta_x);
	} else {
		slope_minus.rho = ((grid[cellnum].state.rho - state_minus.rho)/delta_y); // was 2*delta_y, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((grid[cellnum].state.rhou - state_minus.rhou)/delta_y);
		slope_minus.rhov = ((grid[cellnum].state.rhov - state_minus.rhov)/delta_y);
		slope_minus.E = ((grid[cellnum].state.E - state_minus.E)/delta_y);
		slope_plus.rho = ((state_plus.rho - grid[cellnum].state.rho)/delta_y);
		slope_plus.rhou = ((state_plus.rhou - grid[cellnum].state.rhou)/delta_y);
		slope_plus.rhov = ((state_plus.rhov - grid[cellnum].state.rhov)/delta_y);
		slope_plus.E = ((state_plus.E - grid[cellnum].state.E)/delta_y);
	}
}


void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number) {
	switch(limiter_number) {
		case 1: // Harmonic
			harmonic_limiter(limiter_value, slope_minus, slope_plus);
			break;
//		default: do default stuff
//			break;
	}
}


void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus) {
	double delta = pow(10,-12);
	limiter_value.rho = (std::abs(slope_minus.rho)*slope_plus.rho + slope_minus.rho*std::abs(slope_plus.rho))/(std::abs(slope_minus.rho) + std::abs(slope_plus.rho) + delta);
	limiter_value.rhou = (std::abs(slope_minus.rhou)*slope_plus.rhou + slope_minus.rhou*std::abs(slope_plus.rhou))/(std::abs(slope_minus.rhou) + std::abs(slope_plus.rhou) + delta);
	limiter_value.rhov = (std::abs(slope_minus.rhov)*slope_plus.rhov + slope_minus.rhou*std::abs(slope_plus.rhov))/(std::abs(slope_minus.rhov) + std::abs(slope_plus.rhov) + delta);
	limiter_value.E = (std::abs(slope_minus.E)*slope_plus.E + slope_minus.E*std::abs(slope_plus.E))/(std::abs(slope_minus.E) + std::abs(slope_plus.E) + delta);
}


void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, int direction, double gamma) {
// flux is for the (i-1/2 and i+1/2) faces [or (j-1/2 and j+1/2) faces]
	TDstate center_flux; //computed in line below
	compute_center_flux(center_flux, center_state, direction, gamma);
//	std::cout << "Non-limited flux: " << center_flux.rho << '\n';

	assert(flux.size() == 2);

//	left edge rho flux = (center rho flux)/(center rho)* [(center_rho) +/- 0.5*(plus half rho slope)*(rho limiter)*(1 - local CFL number)]
	flux[0].rho = ((center_flux.rho/center_state.rho)*(center_state.rho - 0.5*(slope_plus.rho*limiter_value.rho)));
	flux[0].rhou = ((center_flux.rhou/center_state.rhou)*(center_state.rhou - 0.5*(slope_plus.rhou*limiter_value.rhou)));
	flux[0].rhov = ((center_flux.rhov/center_state.rhov)*(center_state.rhov - 0.5*(slope_plus.rhov*limiter_value.rhov)));
	flux[0].E = ((center_flux.E/center_state.E)*(center_state.E - 0.5*(slope_plus.E*limiter_value.E)));

	flux[1].rho = ((center_flux.rho/center_state.rho)*(center_state.rho + 0.5*(slope_plus.rho*limiter_value.rho)));
	flux[1].rhou = ((center_flux.rhou/center_state.rhou)*(center_state.rhou + 0.5*(slope_plus.rhou*limiter_value.rhou)));
	flux[1].rhov = ((center_flux.rhov/center_state.rhov)*(center_state.rhov + 0.5*(slope_plus.rhov*limiter_value.rhov)));
	flux[1].E = ((center_flux.E/center_state.E)*(center_state.E + 0.5*(slope_plus.E*limiter_value.E)));
}

void compute_center_flux(TDstate &center_flux, TDstate center_state, int direction, double gamma) {
	double U = center_state.rhou/center_state.rho;
	double V = center_state.rhov/center_state.rho;
	double pressure = compute_pressure_2D(center_state, gamma);
	
	if (direction == 0) {
		center_flux.rho = center_state.rhou;
		center_flux.rhou = center_state.rho*pow(U,2) + pressure;
		center_flux.rhov = center_state.rho*U*V;
		center_flux.E = U*(center_state.E + pressure);
	} else { //direction == 1, bottom/top
		center_flux.rho = center_state.rhou;
		center_flux.rhou = center_state.rho*U*V;		
		center_flux.rhov = center_state.rho*pow(V,2) + pressure;
		center_flux.E = V*(center_state.E + pressure);
	}
}

void compute_halfway_state(unsigned int& cellnum, TDstate& Uph, const TDstate& current_U, const std::vector<double>& F, const std::vector<double>& G, double dt, double dx, double dy) {
	Uph[cellnum].rho = (current_U.rho - dt/(2*dx)*(F[1].rho-F[0].rho) - dt/(2*dy)*(G[1].rho-G[0].rho));
	Uph[cellnum].rhou = (current_U.rhou - dt/(2*dx)*(F[1].rhou-F[0].rhou) - dt/(2*dy)*(G[1].rhou-G[0].rhou));
	Uph[cellnum].rhov = (current_U.rhov - dt/(2*dx)*(F[1].rhov-F[0].rhov) - dt/(2*dy)*(G[1].rhov-G[0].rhov));
	Uph[cellnum].E = (current_U.E - dt/(2*dx)*(F[1].E-F[0].E) - dt/(2*dy)*(G[1].E-G[0].E));
}

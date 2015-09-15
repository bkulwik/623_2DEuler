#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include "file_header.h"
#include "helper_functions.h"

void FluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellnum, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double CFL, double& max_wavespeed, std::vector<TDstate>& U, bool debug);

// Define helper functions!
bool isedge(const std::vector<cell>& grid, int cellnum);

// void compute_flux(std::vector<double> &fluxph, std::vector<double> &fluxmh, std::vector<double> slope, TDstate state, double delta, std::string direction, double gamma);

// This is not math-y, just a formality of grabbing the left/right (state_minus) or bottom/top (state_plus) states from the grid input.
void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const std::vector<cell>& grid, int cellnum, int direction);

// Computes the "geometric" slope vector for all conserved quantities between the bottom/left and top/right cells
void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y, TDstate& U_center);

// This computes the double-value "limiter value" based on the bottom/left and top/right states for a given cell. The limiter value is then used to calculate the flux
void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number);

// Computes the harmonic limiter given r, the ratio of minus to plus slopes
void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus);

void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, const TDstate& slope_minus, int direction, double gamma, double CFL);

// Computes the value of flux in the center cell, not along the edges but in the cell
void compute_cell_flux(TDstate &flux, const TDstate& state, int direction, double gamma);

// Reconstructs the solution from t to t+dt/2
void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double dt, double dx, double dy);

int main() {

// Define Variables
	std::string parameter_filename = "solver_parms.txt";
	std::string input_filename = "test4.bkcfd";
	std::string output_filename = "test4_output.txt";
	std::vector<double> parameters;
	std::vector<double> Fiph, Fimh, Giph, Gimh;
	double CFL, dt, delta_x, delta_y, min_delta_x, min_delta_y, min_delta;
	TDstate state_minus, state_plus, slope_minus, slope_plus, limiter_value;
	TDstate state_leftcell, state_rightcell, state_bottomcell, state_topcell;
	std::vector<TDstate> F (2);
	std::vector<TDstate> G (2);

	double thresh = pow(10,-8); // 1E-8 tolerance for riemann solver - velocity
	double gamma = 1.4;
	bool debug = false;
	bool save_timesteps = true;	

// 1 = harmonic mean for slopes, 2 = first order upwind (no limiter) ...
	int limiter_number = 2; 

	std::vector<cell> grid;

	// Read parameters from input file defined by filename string
	parameters = read_parameters(parameter_filename);
	CFL = parameters.at(0);	
	int num_timesteps = parameters.at(1);
	

// Read initial file defined from Matlab with cell edges and connections
// Read the initial condition file from MATLAB with [rho rho*u rho*v e] defined for each cell
	read_grid(input_filename, grid, gamma);
	
	min_delta_x = gridmin(grid, 'x');
	min_delta_y = gridmin(grid, 'y');
	min_delta = std::min(min_delta_x,min_delta_y);
/*
	for(auto input: grid) {
		std::cout << input.cellnumber << '\n';
		std::cout << "cornerlocs_x: " << input.cornerlocs_x[0] << ' ';
		std::cout << input.cornerlocs_x[1] << ' ';
		std::cout << input.cornerlocs_x[2] << ' ';
		std::cout << input.cornerlocs_x[3] << '\n';
		std::cout << "cornerlocs_y: " << input.cornerlocs_y[0] << ' ';
		std::cout << input.cornerlocs_y[1] << ' ';
		std::cout << input.cornerlocs_y[2] << ' ';
		std::cout << input.cornerlocs_y[3] << '\n';
		std::cout << "adjacent_cells: " << input.adjacent_cells[0] << ' ';
		std::cout << input.adjacent_cells[1] << ' ';
		std::cout << input.adjacent_cells[2] << ' ';
		std::cout << input.adjacent_cells[3] << '\n';
		std::cout << "is_edgecell: " << input.is_edgecell << '\n';
		std::cout << "state: " << input.state.rho << ' ';
		std::cout << input.state.rhou << ' ';
		std::cout << input.state.rhov << ' ';
		std::cout << input.state.E << '\n' << '\n';
	} 
*/
	std::vector<TDstate> U;
	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) {
		U.push_back(grid[cellnum].state);
	}

	double max_wavespeed = 0; 	
	for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
		if (grid[cellnum].is_edgecell == 0) { // Cell is not a ghost cell/edge cell
         FluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellnum, limiter_value, F, G, CFL, max_wavespeed, U, debug);
		}
	}
	dt = CFL*min_delta/max_wavespeed;

// Start calculation in for loop going to final time
	TDstate Uph;
	std::vector<TDstate> Up1 (grid.size());
	for (int timestep = 0; timestep < num_timesteps; timestep++) { // for each timestep
		std::cout << "Timestep " << timestep << " of " << num_timesteps << '\n';
		if (timestep > 0) {
			U = Up1;
		}
		if (debug) {
			outputU(U);
		}

// Go through each cell, calculate fluxes, update to new piecewise linear state
		for (unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) { // For each cell
			if (grid[cellnum].is_edgecell == 0)  { // Cell is not a ghost cell/edge cell
   	   	FluxSolver(grid, delta_x, delta_y, state_minus, state_plus, gamma, slope_plus, slope_minus, limiter_number, cellnum, limiter_value, F, G, CFL, max_wavespeed, U, debug);
			// Now, use the fluxes calculated above to assign a new state, Uph
				compute_halfway_state(Uph, U[cellnum], F, G, dt, delta_x, delta_y);
		
				if (debug) {				
				std::cout << "Cell number: " << cellnum << '\n';
				std::cout << "bottom/top Limiter - " << limiter_value.rho << " " << limiter_value.rhou << " " << limiter_value.rhov << " "  << limiter_value.E << '\n';
				std::cout << "Current State: " << U[cellnum].rho << " " << U[cellnum].rhou << " " << U[cellnum].rhov << " " << U[cellnum].E << '\n';
//				std::cout << "F, cell number " << cellnum << ": " << '\n';
//				std::cout << F[0].rho << " " << F[0].rhou << " " << F[0].rhov << " " << F[0].E << '\n';
//				std::cout << F[1].rho << " " << F[1].rhou << " " << F[1].rhov << " " << F[1].E << '\n';
	
//				std::cout << "G, cell number: " << cellnum << '\n';
//				std::cout << G[0].rho << " " << G[0].rhou << " " << G[0].rhov << " " << G[0].E << '\n';
//				std::cout << G[1].rho << " " << G[1].rhou << " " << G[1].rhov << " " << G[1].E << '\n';
	
				std::cout << "Cell " << cellnum << " half-state: " <<  Uph.rho << " " << Uph.rhou << " " << Uph.rhov << " " << Uph.E << '\n';
				}

// Solve riemann problem on edge for euler flux
				TDstate state_center_2D = U[cellnum];
				TDstate tempflux;
				ODstate state_leftcell_1D, center_lr, center_bt, state_rightcell_1D, state_bottomcell_1D, state_topcell_1D;
				
				assign_edge_state(state_leftcell, state_rightcell, U, grid, cellnum, 0); // The left and right cell TDstates of the current cell
				assign_edge_state(state_bottomcell, state_topcell, U, grid, cellnum, 1); // The bottom and top cell TDstates of the current cell
				
// All this mess is because the riemann problem is 1D and we store 2D states.
// So we have to make new, temporary 1D states and fixes (all the doubles) to make the RP work.
				center_lr.rho = U[cellnum].rho;
				center_lr.rhou = U[cellnum].rhou;
				double center_lr_parallelvel = U[cellnum].rhov/U[cellnum].rho; 
				center_lr.E = U[cellnum].E;

				center_bt.rho = U[cellnum].rho;
				double center_bt_parallelvel = U[cellnum].rhou/U[cellnum].rho; 
				center_bt.rhou = U[cellnum].rhov;	
				center_bt.E = U[cellnum].E;

				state_leftcell_1D.rho = state_leftcell.rho;
				state_leftcell_1D.rhou = state_leftcell.rhou;
				double leftcell_parallelvel = state_leftcell.rhov/state_leftcell.rho;
				state_leftcell_1D.E = state_leftcell.E;

				state_rightcell_1D.rho = state_rightcell.rho; 
				state_rightcell_1D.rhou = state_rightcell.rhou;
				double rightcell_parallelvel = state_rightcell.rhov/state_rightcell.rho; 
				state_rightcell_1D.E = state_rightcell.E;

				state_bottomcell_1D.rho = state_bottomcell.rho;
				double bottomcell_parallelvel = state_bottomcell.rhou/state_bottomcell.rho;
				state_bottomcell_1D.rhou = state_bottomcell.rhov;
				state_bottomcell_1D.E = state_bottomcell.E;

				state_topcell_1D.rho = state_topcell.rho;
				double topcell_parallelvel = state_topcell.rhou/state_topcell.rho;
				state_topcell_1D.rhou = state_topcell.rhov;
				state_topcell_1D.E = state_topcell.E;


// Need to read in overall velocity magnitude to Exact Riemann Solver so we can get an accurate pressure reading!
// Reconstruct, from the state_[  ] ODstates, the TDstates corresponding assuming that the non-OD direction velocity is constant
// Then compute the flux from these new TDstates (4x), and update solution from t+1/2 to t+1 !!!
				ODstate state_leftedge_1D = Exact_Riemann_Solver(state_leftcell_1D, center_lr, leftcell_parallelvel, center_lr_parallelvel, thresh, gamma, debug);
				ODstate state_rightedge_1D = Exact_Riemann_Solver(center_lr, state_rightcell_1D, center_lr_parallelvel, rightcell_parallelvel, thresh, gamma, debug);
				ODstate state_bottomedge_1D = Exact_Riemann_Solver(state_bottomcell_1D, center_bt, bottomcell_parallelvel, center_bt_parallelvel, thresh, gamma, debug);
				ODstate state_topedge_1D = Exact_Riemann_Solver(center_bt, state_topcell_1D, center_bt_parallelvel, topcell_parallelvel, thresh, gamma, debug);


				TDstate state_leftedge_2D, state_rightedge_2D, state_bottomedge_2D, state_topedge_2D;


				state_leftedge_2D.rho = state_leftedge_1D.rho;
				state_leftedge_2D.rhou = state_leftedge_1D.rhou;
				state_leftedge_2D.rhov = state_center_2D.rhov;
				state_leftedge_2D.E = state_leftedge_1D.E;
				
				state_rightedge_2D.rho = state_rightedge_1D.rho;
				state_rightedge_2D.rhou = state_rightedge_1D.rhou;
				state_rightedge_2D.rhov = state_center_2D.rhov;
				state_rightedge_2D.E = state_rightedge_1D.E;

				state_bottomedge_2D.rho = state_bottomedge_1D.rho;
				state_bottomedge_2D.rhou = state_center_2D.rhou; 
				state_bottomedge_2D.rhov = state_bottomedge_1D.rhou;
				state_bottomedge_2D.E = state_bottomedge_1D.E;

				state_topedge_2D.rho = state_topedge_1D.rho;
				state_topedge_2D.rhou = state_center_2D.rhou; 
				state_topedge_2D.rhov = state_topedge_1D.rhou;
				state_topedge_2D.E = state_topedge_1D.E;

				compute_cell_flux(tempflux, state_leftedge_2D, 0, gamma);
				F[0] = tempflux;
				compute_cell_flux(tempflux, state_rightedge_2D, 0, gamma);
				F[1] = tempflux;
				compute_cell_flux(tempflux, state_bottomedge_2D, 1, gamma);
				G[0] = tempflux;
				compute_cell_flux(tempflux, state_topedge_2D, 1, gamma);
				G[1] = tempflux;
//Use all these fluxes to define updated state on each cell 
				compute_halfway_state(Up1[cellnum], Uph, F, G, dt, delta_x, delta_y); 
	
				if (debug) {
					std::cout << "Cell " << cellnum << " Updated State : " <<  Up1[cellnum].rho << " " << Up1[cellnum].rhou << " " << Up1[cellnum].rhov << " " << Up1[cellnum].E << '\n' << '\n';
				}
			} // end of if (!isedge(grid,cell)) 
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) (first one)

		// Now need to redefine the ghost cells to the correct values
		for(unsigned int cellnum = 0; cellnum < grid.size(); ++cellnum) {
			if (grid[cellnum].is_edgecell == 1) { // It is a ghost cell 
				int special_value = -1;		
				int num_edges = 0;		
				std::vector<int> edge_locations;
				for(int adjacentcell = 0; adjacentcell < 4; ++adjacentcell) {
					if (grid[cellnum].adjacent_cells[adjacentcell] == special_value) {
						edge_locations.push_back(adjacentcell);
						num_edges++;
					}
				}
		//		if (num_edges >= 2) {
	//				continue; // Sorry bro, we don't care about updating you!
//				} else if (num_edges == 1) {
					if (edge_locations[0] < 2) {
						edge_locations[0] = (edge_locations[0] + 2);
					} else {
						edge_locations[0] = (edge_locations[0] - 2);					
					}
				
					Up1[cellnum] = Up1[grid[cellnum].adjacent_cells[edge_locations[0]]];
					if (edge_locations[0] >= 2) {
						edge_locations[0] = (edge_locations[0] - 2);
					}
					assert((edge_locations[0] == 0) || (edge_locations[0] == 1));
					if (edge_locations[0] == 0) { // left-right
						if (Up1[cellnum].rhou == 0) {
							continue;
						} else {
							Up1[cellnum].rhou = -(Up1[cellnum].rhou);
						}
					} else if (edge_locations[0] == 1) { //up-down
						if (Up1[cellnum].rhou == 0) {
							continue;
						} else {
							Up1[cellnum].rhov = -(Up1[cellnum].rhov);
						}
					}
				//}
			} // end of if (isedge(grid,cell)) 
		} // end of for(unsigned int cell = 0; cell < grid.size(); ++cell) (second one)
		if (save_timesteps) {
			write_to_file(grid,Up1,output_filename);
		}
	} // end of for(double t = 0; t <= 2*dt; t++) 
return 0;
}



// This is the land of the helper functions! ---------------------------------------------------

void FluxSolver(const std::vector<cell>& grid, double& delta_x, double& delta_y, TDstate& state_minus, TDstate& state_plus, double gamma, TDstate& slope_plus, TDstate& slope_minus, int limiter_number, int cellnum, TDstate& limiter_value, std::vector<TDstate>& F, std::vector<TDstate>& G, double CFL, double& max_wavespeed, std::vector<TDstate>& U, bool debug) {

	delta_x = (vectormax(grid[cellnum].cornerlocs_x) - vectormin(grid[cellnum].cornerlocs_x));
	delta_y = (vectormax(grid[cellnum].cornerlocs_y) - vectormin(grid[cellnum].cornerlocs_y));

    // Calculate the limited flux of each conserved quantity for initial t -> t+1/2 update
	for (int direction = 0; direction < 2; ++direction) { // For the left/right direction and up/down direction
		// Assign the bottom/left and the top/right states for each cell
		assign_edge_state(state_minus, state_plus, U, grid, cellnum, direction);	
						
		// Calculate the flux of each conserved variable, in each direction, on each face, of that cell.
		// Flux needs to be calculated by first finding the limiter (direction-dependent), then using that limiter to calculate the directional flux, then using the two directional fluxes to update from u_t to u_t+1/2	
			
		compute_slope(slope_minus, slope_plus, direction, state_minus, state_plus, delta_x, delta_y, U[cellnum]);
		if ((direction == 1) && (debug)) {
			std::cout << "Slope Minus: " << slope_minus.rho << " " << slope_minus.rhou << " " << slope_minus.rhov << " " << slope_minus.E << '\n';
			std::cout << "Slope Plus: " << slope_plus.rho << " " << slope_plus.rhou << " " << slope_plus.rhov << " " << slope_plus.E << '\n';
		}
		compute_limiter(limiter_value, slope_minus, slope_plus, limiter_number);

		if (direction == 0) { // Compute F, x-fluxes
			compute_limited_flux(F, limiter_value, U[cellnum], slope_plus, slope_minus, direction, gamma, CFL);
//			std::cout << "rho limiter x: " << limiter_value.rho << '\n';
		} else { // Compute G, y-fluxes
//			std::cout << "rho limiter y: " << limiter_value.rho << '\n';
			compute_limited_flux(G, limiter_value, U[cellnum], slope_plus, slope_minus, direction, gamma, CFL);
		}
	} // end of for (int direction = 0; direction < 2; ++direction)

	max_wavespeed_calculator_riemann(max_wavespeed, U[cellnum], gamma);
}


bool isedge(const std::vector<cell>& grid, int cellnum) {
	int special_value = -1;
	if ((grid[cellnum].adjacent_cells[0] == special_value) || (grid[cellnum].adjacent_cells[1] == special_value) || (grid[cellnum].adjacent_cells[2] == special_value)  || (grid[cellnum].adjacent_cells[3] == special_value)) {
		return(1);
	} else {
		return(0);
	}
}


void assign_edge_state(TDstate &state_minus, TDstate &state_plus, const std::vector<TDstate>& U, const std::vector<cell>& grid, int cellnum, int direction) {
	int adjacent_bl, adjacent_tr;

	adjacent_bl = (grid[cellnum].adjacent_cells[direction]);
	adjacent_tr = (grid[cellnum].adjacent_cells[direction+2]);

	state_minus = U[adjacent_bl];
	state_plus = U[adjacent_tr];
}


void compute_slope(TDstate &slope_minus, TDstate &slope_plus, int direction, TDstate state_minus, TDstate state_plus, double delta_x, double delta_y, TDstate& U_center) {
	if (direction == 0) { //left-right
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta_x); // was 2*delta_x, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta_x);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta_x);
		slope_minus.E = ((U_center.E - state_minus.E)/delta_x);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta_x);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta_x);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta_x);
		slope_plus.E = ((state_plus.E - U_center.E)/delta_x);
	} else {
		slope_minus.rho = ((U_center.rho - state_minus.rho)/delta_y); // was 2*delta_y, changed to just delta_x from 1/21/15 notes pg.2
		slope_minus.rhou = ((U_center.rhou - state_minus.rhou)/delta_y);
		slope_minus.rhov = ((U_center.rhov - state_minus.rhov)/delta_y);
		slope_minus.E = ((U_center.E - state_minus.E)/delta_y);

		slope_plus.rho = ((state_plus.rho - U_center.rho)/delta_y);
		slope_plus.rhou = ((state_plus.rhou - U_center.rhou)/delta_y);
		slope_plus.rhov = ((state_plus.rhov - U_center.rhov)/delta_y);
		slope_plus.E = ((state_plus.E - U_center.E)/delta_y);
	}
}


void compute_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus, int limiter_number) {
	switch(limiter_number) {
		case 1: // Harmonic
			harmonic_limiter(limiter_value, slope_minus, slope_plus);
			break;
		case 2: // First order upwind
			limiter_value.rho = 0;
			limiter_value.rhou = 0;
			limiter_value.rhov = 0;
			limiter_value.E = 0;
//		default: do default stuff
//			break;
	}
}


void harmonic_limiter(TDstate &limiter_value, TDstate slope_minus, TDstate slope_plus) {
	double delta = pow(10,-12);
	limiter_value.rho = (std::abs(slope_minus.rho)*slope_plus.rho + slope_minus.rho*std::abs(slope_plus.rho))/(std::abs(slope_minus.rho) + std::abs(slope_plus.rho) + delta);
	limiter_value.rhou = (std::abs(slope_minus.rhou)*slope_plus.rhou + slope_minus.rhou*std::abs(slope_plus.rhou))/(std::abs(slope_minus.rhou) + std::abs(slope_plus.rhou) + delta);
	limiter_value.rhov = (std::abs(slope_minus.rhov)*slope_plus.rhov + slope_minus.rhov*std::abs(slope_plus.rhov))/(std::abs(slope_minus.rhov) + std::abs(slope_plus.rhov) + delta);
	limiter_value.E = (std::abs(slope_minus.E)*slope_plus.E + slope_minus.E*std::abs(slope_plus.E))/(std::abs(slope_minus.E) + std::abs(slope_plus.E) + delta);
}

void compute_limited_flux(std::vector<TDstate>& flux, const TDstate& limiter_value, const TDstate& center_state, const TDstate& slope_plus, const TDstate& slope_minus, int direction, double gamma, double CFL) {
// flux is for the (i-1/2 and i+1/2) faces [or (j-1/2 and j+1/2) faces]
	TDstate minus_limited_state, plus_limited_state, tempflux;

	assert(flux.size() == 2);

	minus_limited_state.rho = center_state.rho - 0.5*slope_minus.rho*limiter_value.rho*(1-CFL);
	minus_limited_state.rhou = center_state.rhou - 0.5*slope_minus.rhou*limiter_value.rhou*(1-CFL);
	minus_limited_state.rhov = center_state.rhov - 0.5*slope_minus.rhov*limiter_value.rhov*(1-CFL);
	minus_limited_state.E = center_state.E - 0.5*slope_minus.E*limiter_value.E*(1-CFL);

	plus_limited_state.rho = center_state.rho + 0.5*slope_plus.rho*limiter_value.rho*(1-CFL);
	plus_limited_state.rhou = center_state.rhou + 0.5*slope_plus.rhou*limiter_value.rhou*(1-CFL);
	plus_limited_state.rhov = center_state.rhov + 0.5*slope_plus.rhov*limiter_value.rhov*(1-CFL);
	plus_limited_state.E = center_state.E + 0.5*slope_plus.E*limiter_value.E*(1-CFL);

	compute_cell_flux(tempflux, minus_limited_state, direction, gamma);
	flux[0] = (tempflux);
	compute_cell_flux(tempflux, plus_limited_state, direction, gamma);
	flux[1] = (tempflux);
}

void compute_cell_flux(TDstate &flux, const TDstate& state, int direction, double gamma) {
	assert(state.rho > 0);
	assert(state.E > 0);

	double U = state.rhou/state.rho;
	double V = state.rhov/state.rho;
	double pressure = compute_pressure_2D(state, gamma);
	
	if (direction == 0) {
		flux.rho = state.rhou;
		flux.rhou = state.rho*pow(U,2) + pressure;
		if ((U == 0) || (V == 0)) {
			flux.rhov = 0;
		} else {
			flux.rhov = state.rho*U*V;
		}
		flux.E = U*(state.E + pressure);
	} else { //direction == 1, bottom/top
		flux.rho = state.rhov;
		if ((U == 0) || (V == 0)) {
			flux.rhou = 0;
		} else {
			flux.rhou = state.rho*U*V;
		}
//		std::cout << "In CCF: rho,V, V^2, P = " << state.rho << " " << V << " " << pow(V,2) << " " << pressure << '\n';
		flux.rhov = state.rho*pow(V,2) + pressure;
		flux.E = V*(state.E + pressure);
	}
}

void compute_halfway_state(TDstate& Uph, const TDstate& current_U, const std::vector<TDstate>& F, const std::vector<TDstate>& G, double dt, double dx, double dy) {
	Uph.rho = (current_U.rho - dt/(2*dx)*(F[1].rho-F[0].rho) - dt/(2*dy)*(G[1].rho-G[0].rho));
	Uph.rhou = (current_U.rhou - dt/(2*dx)*(F[1].rhou-F[0].rhou) - dt/(2*dy)*(G[1].rhou-G[0].rhou));
	Uph.rhov = (current_U.rhov - dt/(2*dx)*(F[1].rhov-F[0].rhov) - dt/(2*dy)*(G[1].rhov-G[0].rhov));
	Uph.E = (current_U.E - dt/(2*dx)*(F[1].E-F[0].E) - dt/(2*dy)*(G[1].E-G[0].E));
}

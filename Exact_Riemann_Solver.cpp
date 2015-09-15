#include "file_header.h"
#include "helper_functions.h"
#include <iostream>
#include <cmath>

//left is the [rho rho*u E] for the 'left' state - the cell on the bottom or left
//right_state is the [rho rho*u E] for the 'right' state - the cell on the top or right
//
// This entire function was made using the guidance and notation of:
// http://www.ifm.umich.mx/guzman/Grupo/Educacion_4.pdf

struct wavetype {
	int left;
	int right;
} type;

ODstate compute_zerostate(const std::vector<double>& wavespeeds, const ODstate& left, const ODstate& U_star_l, const ODstate& U_star_r, const ODstate& right, double gamma, double a_l, double a_r, double left_parallelvel, double right_parallelvel);

ODstate Exact_Riemann_Solver(ODstate left, ODstate right, double left_parallelvel, double right_parallelvel, double thresh, double gamma, bool debug) {

	std::cout.precision(15);

	double v_l = left.rhou/left.rho;
	double v_r = right.rhou/right.rho;
// Use compute_pressure_2D to find the actual pressure - with 2D velocity
//----------------------------------------------------------------
	TDstate left_TD, right_TD;

	left_TD.rho = left.rho;
	left_TD.rhou = left.rhou;
	left_TD.rhov = left.rho*left_parallelvel;
	left_TD.E = left.E;	

	right_TD.rho = right.rho;
	right_TD.rhou = right.rhou;
	right_TD.rhov = right.rho*right_parallelvel;
	right_TD.E = right.E;

	double p_l = compute_pressure_2D(left_TD, gamma);
	double p_r = compute_pressure_2D(right_TD, gamma);
//----------------------------------------------------------------

	double a_l = sqrt(gamma*p_l/left.rho);
	double a_r = sqrt(gamma*p_r/right.rho);
	double a_star_l, a_star_r;

	std::vector<double> wavespeeds_l, wavespeeds_r, wavespeeds_star, wavespeeds;

	double delta, p_guess;
	int iter = 1;
	double limiter = 1;
	std::vector<double> p_guess_temp(3), rho_l_star (3), rho_r_star (3), v_l_star (3), v_r_star (3);
	std::vector<double> E_l_star (3), E_r_star (3), p_l_star (3), p_r_star (3);
	std::vector<double> error (3);
	double derrordp, A, B;
	
	ODstate left_new, right_new, U_star_l, U_star_r;

	//Determine which types to solve (rarefaction(1), shock(2), or nothing)
	//We are allowed a density discontinuity but not a pressure or velocity discontinuity
	if (p_l > p_r) {
		type.left = 1;
		type.right = 2;
	}
	else if (p_l < p_r) {
		type.left = 2;
		type.right = 1;
	}
	else if (v_l < v_r) {
		type.left = 1;
		type.right = 1;
	}
	else if (v_l > v_r) {
		type.left = 2;
		type.right = 2;
	}
	else {
		type.left = 0;
		type.right = 0;
		left_new = left;
		right_new = right;
		// Not sure how to handle this - just check to see if we maintain constant state I guess...
		if (v_l < 0) {
			return (right);
		}
		else {
			return(left);
		}
	}

	delta = std::max(p_l, p_r)*0.000001; //Used to find the derivative of the error with respect to p, in order to iterate on p in the star state effectively and find the actual solution

	p_guess = 0.5*(p_l + p_r); //Guess for star state pressure

	while(1) {
//		std::cout << '\n' << "Iteration Number " << iter << '\n';
		
//		std::cout << "Pressure Guess: " << p_guess << '\n';
		p_guess_temp.at(0) = (p_guess - delta);
		p_guess_temp.at(1) = p_guess;
		p_guess_temp.at(2) = (p_guess + delta);
		
/*		std::cout << "p_guess_temp " << '\n';
		for(auto iF:p_guess_temp) {
				std::cout << iF << '\n';
		}
*/
		
		if (type.left == 1) { //The left state undergoes a rarefaction
			for (int it = 0; it < 3; ++it) {
				v_l_star[it] = (v_l - 2*a_l/(gamma-1)*(pow((p_guess_temp[it]/p_l), ((gamma-1)/(2*gamma))) - 1));
				rho_l_star[it] = (left.rho*pow((p_guess_temp[it]/p_l),(1/gamma)));
//				rho_l_star[it] = left.rho*pow((2/(gamma+1) + (gamma-1)/(a_l*(gamma+1))*(left.rhou/left.rho)),(2/(gamma-1)));
				E_l_star[it] = ((p_guess_temp[it]/(gamma-1)) - rho_l_star[it]*pow(v_l_star[it],2)/2);			
			}
		}
		else { //The left state undergoes a shock
			A = 2/((gamma-1)*left.rho);
			B = (gamma - 1)/(gamma + 1)*p_l;
			for (int it = 0; it < 3; ++it) {
				v_l_star[it] = (v_l - (p_guess_temp[it]-p_l)*sqrt(A/(p_guess_temp[it] + B)));
				rho_l_star[it] = (left.rho*((p_l*(gamma-1) + p_guess_temp[it]*(gamma+1))/(p_guess_temp[it]*(gamma-1) + p_l*(gamma+1))));
				E_l_star[it] = ((p_guess_temp[it]/(gamma-1)) - rho_l_star[it]*pow(v_l_star[it],2)/2);
			}
		}

		if (type.right == 1) { //The right state undergoes a rarefaction
			for (int it = 0; it < 3; ++it) {
				v_r_star[it] = (v_r - 2*a_r/(gamma-1)*(1 - pow((p_guess_temp[it]/p_r),((gamma-1)/(2*gamma)))));
				rho_r_star[it] = (right.rho*pow((p_guess_temp[it]/p_r),(1/gamma)));
//				rho_r_star[it] = right.rho*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),(2/(gamma-1))); 
				E_r_star[it] = ((p_guess_temp[it]/(gamma-1)) - rho_r_star[it]*pow(v_r_star[it],2)/2);			
			}
		}
		else { //The right state undergoes a shock
			A = 2/((gamma-1)*right.rho);
			B = (gamma - 1)/(gamma + 1)*p_r;
			for (int it = 0; it < 3; ++it) {
				v_r_star[it] = (v_r + (p_guess_temp[it]-p_r)*sqrt(A/(p_guess_temp[it] + B)));
				rho_r_star[it] = (right.rho*((p_r*(gamma-1) + p_guess_temp[it]*(gamma+1))/(p_guess_temp[it]*(gamma-1) + p_r*(gamma+1))));
				E_r_star[it] = ((p_guess_temp[it]/(gamma-1)) - rho_r_star[it]*pow(v_r_star[it],2)/2);
			}
		}

 /* Error Checking
		std::cout << "v_r_star " << '\n';
		for(auto iF:v_r_star) {
				std::cout << iF << '\n';
		}
		std::cout << "v_l_star " << '\n';
		for(auto iF:v_l_star) {
				std::cout << iF << '\n';
		}
		std::cout << "rho_r_star " << '\n';
		for(auto iF:rho_r_star) {
				std::cout << iF << '\n';
		}
		std::cout << "rho_l_star " << '\n';
		for(auto iF:rho_l_star) {
				std::cout << iF << '\n';
		}
		std::cout << "E_r_star " << '\n';
		for(auto iF:E_r_star) {
				std::cout << iF << '\n';
		}
		std::cout << "E_l_star " << '\n';
		for(auto iF:E_l_star) {
				std::cout << iF << '\n';
		}
 End */ 


// Start using the Newton's method solver here
		for (int it = 0; it < 3; ++it) {
			error[it] = (v_l_star[it] - v_r_star[it]);
		}

		//std::cout << "error = " << error.at(1) << '\n';

		//Rate of error change with p
		derrordp = (error.at(2) - error.at(0))/(2*delta);
		if (std::abs(error.at(1)) > thresh) {
/* Error Checking
			std::cout << "Limiter: " << limiter << '\n';
			std::cout << "derrordp: " << derrordp << '\n';
			std::cout << "Pressure Guess change: " << limiter*error.at(1)/derrordp << '\n';
End */
			iter++;
			p_guess = p_guess - limiter*error.at(1)/derrordp;
			if (p_guess <= 0) {
				p_guess = delta;
			}
		}
		else { //Yey, the calculation converged on a pressure for the star state! And therefore a velocity, and therefore an energy!
			U_star_l.rho = rho_l_star[1]; 
			U_star_l.rhou = (U_star_l.rho*v_l_star[1]);
			U_star_l.E = (p_guess_temp[1]/(gamma-1)-(pow(U_star_l.rhou,2)/(U_star_l.rho*2))); 
			U_star_l.pressure = p_guess_temp[1];
			a_star_l = sqrt(gamma*p_guess_temp[1]/U_star_l.rho); //sqrt(gamma*P/rho)

			U_star_r.rho = rho_r_star[1]; 
			U_star_r.rhou = (U_star_r.rho*v_r_star[1]);
			U_star_r.E = (p_guess_temp[1]/(gamma-1)-(pow(U_star_r.rhou,2)/(U_star_r.rho*2))); 
			U_star_r.pressure = p_guess_temp[1];
			a_star_r = sqrt(gamma*p_guess_temp[1]/U_star_r.rho); //sqrt(gamma*P/rho)

 // Error Checking
	/*		std::cout << '\n' << '\n' << "Yey, we converged!! U Star Left is: " << '\n';
			for(auto ustar:U_star_l) {
				std::cout << ustar << '\n';
			}
			std::cout << "V_l = " << U_star_l.rhou/U_star_l.rho << '\n';
			std::cout << '\n' << "U Star Right is: " << '\n';
			for(auto ustar:U_star_r) {
				std::cout << ustar << '\n';
			}
			std::cout << "V_r = " << U_star_r.rhou/U_star_r.rho << '\n';
*/
//End 

			break;
				// The long term is E
				// E = p/(gamma-1) - rho*u^2/2;W
		}

/*	//	Run checks to make sure solution is converging - this was needed in matlab script but not used here
		if (iter > 50 && num_limiterchanges <= 25) {
			limiter = limiter/2;
			iter = 1;
			num_limiterchanges++;
			std::cout << "Does not converge - decreasing limiter value to " << limiter << '\n';
		}
		else if (iter > 50 && num_limiterchanges > 25) {
			std::cout << "Does not converge, density difference " << rho_r_star.rhou - rho_l_star.rhou << '\n';
		}
*/		
	}

	if (debug) {
	//	std::cout << "Left " ;
	//	outputODstate(left);
	//	std::cout << "Left star " ;
	//	outputODstate(U_star_l);
	//	std::cout << "Right star " ;
	//	outputODstate(U_star_r);
	//	std::cout << "Right " ;	
	//	outputODstate(right);
	}

	// Now compute which state is on the x=0 line

	// First, find and order wavespeeds
	// Note - according to paper, need to use modified wavespeeds!(eq.42) - try with these basic here to see if good solutions
	if (type.left == 1) { //The left state undergoes a rarefaction
		wavespeeds.push_back(left.rhou/left.rho - a_l);
		wavespeeds.push_back(U_star_l.rhou/U_star_l.rho - a_star_l);
	}
	else { //The left state undergoes a shock
		wavespeeds.push_back(left.rhou/left.rho - a_l);
	}

	//average of velocities in left and right star state - should be equal but averaging to reduce numerical error
	wavespeeds.push_back(((U_star_r.rhou/U_star_r.rho) + (U_star_l.rhou/U_star_l.rho))/2); 

	if (type.right == 1) { //The right state undergoes a rarefaction
		wavespeeds.push_back(U_star_r.rhou/U_star_r.rho + a_star_r);
		wavespeeds.push_back(right.rhou/right.rho + a_r);
	}
	else { //The right state undergoes a shock
		wavespeeds.push_back(right.rhou/right.rho + a_r);
	}

	if (debug) {
//		std::cout << "Wavespeeds: ";
//		for(auto input: wavespeeds) {
//			std::cout << input << ", ";
//		}
//		std::cout << '\n';
	}


// Compute and return the state on x = 0
	return (compute_zerostate(wavespeeds, left, U_star_l, U_star_r, right, gamma, a_l, a_r, left_parallelvel, right_parallelvel));	
}




// ------------------------------------------------------------------------------------------

ODstate compute_zerostate(const std::vector<double>& wavespeeds, const ODstate& left, const ODstate& U_star_l, const ODstate& U_star_r, const ODstate& right, double gamma, double a_l, double a_r, double left_parallelvel, double right_parallelvel) {

	int left_wavenumber;
	bool computed = 0;
	ODstate in_raref;
	
	for (unsigned int it = 0; it < wavespeeds.size(); ++it) {
		if (wavespeeds[it] < 0) {
			continue;
		}
		else if (wavespeeds[it] > 0) { // We found the right wave of the center state region 
			computed = 1;
			left_wavenumber = it-1;
		}
	}
	if (!computed) {
		return(right);
	}
	else {
		if (left_wavenumber == 0) { // x=0 to the left of everything
			return (left);
		}
		else if ((left_wavenumber == 1) && type.left == 2) { // left shock and x=0 to right of it but before CS
			return(U_star_l);
		}
		else if ((left_wavenumber == 1) && type.left == 1) { // left rarefaction with x=0 inside it
			in_raref.rho = left.rho*pow((2/(gamma+1) + (gamma-1)/(a_l*(gamma+1))*(left.rhou/left.rho)),(2/(gamma-1)));
			in_raref.rhou = in_raref.rho*(2/(gamma+1)*(a_l + (gamma-1)/2*(left.rhou/left.rho)));

			double press, press_l;
			TDstate press_calculation_state;
			press_calculation_state.rho = left.rho;
			press_calculation_state.rhou = left.rhou; // Order doesn't matter for directionality for this function
			press_calculation_state.rhov = left_parallelvel*left.rho;
			press_calculation_state.E = left.E;

			press_l = compute_pressure_2D(press_calculation_state, gamma);
			press = press_l*pow((2/(gamma+1) + (gamma-1)/(a_l*(gamma+1))*(left.rhou/left.rho)),((2*gamma)/(gamma-1)));


			in_raref.E = (press/(gamma-1) - pow(in_raref.rhou,2)/(2*in_raref.rho));
			return(in_raref);		
		}
		else if ((left_wavenumber == 2) && type.left == 1) { // left rarefaction and x=0 to right of it all but before CS
			return(U_star_l);
		}
		else if ((left_wavenumber == 2) && type.left == 2) { // left shock and x=0 to the right of CS
			return(U_star_r);
		}
		else if ((left_wavenumber == 3) && type.left == 1) { // left rarefaction and x=0 to the right of CS
			return(U_star_r);
		}
		else if ((left_wavenumber == 3) && type.left == 2 && type.right == 1) { // left shock and right rarefaction and x=0 inside right rarefaction
			in_raref.rho = right.rho*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),(2/(gamma-1)));
			in_raref.rhou = in_raref.rho*(2/(gamma+1)*(-a_r + (gamma-1)/2*(right.rhou/right.rho)));

			double press, press_r;
			TDstate press_calculation_state;
			press_calculation_state.rho = right.rho;
			press_calculation_state.rhou = right.rhou; // Order doesn't matter for directionality for this function
			press_calculation_state.rhov = right_parallelvel*right.rho;
			press_calculation_state.E = right.E;

			press_r = compute_pressure_2D(press_calculation_state, gamma);
			press = press_r*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),((2*gamma)/(gamma-1)));

			in_raref.E = (press/(gamma-1) - pow(in_raref.rhou,2)/(2*in_raref.rho));
			return(in_raref);	
		}
//		else if ((left_wavenumber == 3) && type.left == 2 && type.right == 2) { // left shock and right shock and x=0 is to the right of everything - should never be accessed
//			return(right);
//		}
//		else if ((left_wavenumber == 4) && type.left == 1 && type.right == 2) { // left rarefaction and right shock and x=0 is to the right of everything - should never be accessed
//			return(right);
//		}
		else if ((left_wavenumber == 4) && type.left == 1 && type.right == 1) { // left rarefaction and right rarefaction and x=0 is inside right rarefaction
			in_raref.rho = right.rho*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),(2/(gamma-1)));
			in_raref.rhou = in_raref.rho*(2/(gamma+1)*(-a_r + (gamma-1)/2*(right.rhou/right.rho)));

			double press, press_r;
			TDstate press_calculation_state;
			press_calculation_state.rho = right.rho;
			press_calculation_state.rhou = right.rhou; // Order doesn't matter for directionality for this function
			press_calculation_state.rhov = right_parallelvel*right.rho;
			press_calculation_state.E = right.E;

			press_r = compute_pressure_2D(press_calculation_state, gamma);
			press = press_r*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),((2*gamma)/(gamma-1)));

			in_raref.E = (press/(gamma-1) - pow(in_raref.rhou,2)/(2*in_raref.rho));
			return(in_raref);	
		}
//		else if ((left_wavenumber == 4) && type.left == 2) { // left shock and x=0 to the right of everything - should never be accessed!
//			return(right);
//		}
//		else if (left_wavenumber == 5) { // x=0 is to the right of everything - should never be accessed!
//			return(right);
//		}
		else {
			std::cout << '\n' << '\n' << "ERROR: CANNOT FIND WAVE ON X = 0 IN EXACT_RIEMANN_SOLVER.CPP...EXITING" << '\n';
		}
	}
}

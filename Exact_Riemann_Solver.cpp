#include "file_header.h"
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

ODstate compute_zerostate(std::vector<double> wavespeeds, ODstate left, ODstate U_star_l, ODstate U_star_r, ODstate right, double gamma, double a_r, double a_l);

ODstate Exact_Riemann_Solver(ODstate left, ODstate right, double thresh, double gamma) {
	
	std::cout.precision(15);

	double v_l = left.rhou/left.rho;
	double v_r = right.rhou/right.rho;
	double p_l = (gamma-1)*(left.E + pow(left.rhou,2)/(2*left.rho));
	double p_r = (gamma-1)*(right.E + pow(right.rhou,2)/(2*right.rho));
	double a_l = sqrt(gamma*p_l/left.rho);
	double a_r = sqrt(gamma*p_r/right.rho);
	double a_star_l, a_star_r;

	//std::vector<int> type (2);
	std::vector<double> wavespeeds_l, wavespeeds_r, wavespeeds_star, wavespeeds;

	
	double delta, p_guess;
	int iter = 1;
	double limiter = 1;
//	int num_limiterchanges = 0;
	std::vector<double> p_guess_temp(3), rho_l_star (3), rho_r_star (3), v_l_star (3), v_r_star (3);
	std::vector<double> E_l_star (3), E_r_star (3), p_l_star (3), p_r_star (3);
	std::vector<double> error (3);
	double derrordp, A, B;
	
	ODstate left_new, right_new, U_star_l, U_star_r;

	//Determine which types to solve (rarefaction(1), shock(2), or nothing)

	if (left.rho > right.rho) {
		type.left = 1;
		type.right = 2;
	}
	else if (right.rho > left.rho) {
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
		std::cout << '\n' << "Iteration Number " << iter << '\n';
		
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
				v_l_star.at(it) = (v_l - 2*a_l/(gamma-1)*(pow((p_guess_temp.at(it)/p_l), ((gamma-1)/(2*gamma))) - 1));
				rho_l_star.at(it) = (left.rho*pow((p_guess_temp.at(it)/p_l),(1/gamma)));
				E_l_star.at(it) = ((p_guess_temp.at(it)/(gamma-1)) - rho_l_star.at(it)*pow(v_l_star.at(it),2)/2);			
			}
		}
		else { //The left state undergoes a shock
			A = 2/((gamma-1)*left.rho);
			B = (gamma - 1)/(gamma + 1)*p_l;
			for (int it = 0; it < 3; ++it) {
				v_l_star.at(it) = (v_l - (p_guess_temp.at(it)-p_l)*sqrt(A/(p_guess_temp.at(it) + B)));
				rho_l_star.at(it) = (left.rho*((p_l*(gamma-1) + p_guess_temp.at(it)*(gamma+1))/(p_guess_temp.at(it)*(gamma-1) + p_l*(gamma+1))));
				E_l_star.at(it) = ((p_guess_temp.at(it)/(gamma-1)) - rho_l_star.at(it)*pow(v_l_star.at(it),2)/2);
			}
		}

		if (type.right == 1) { //The right state undergoes a rarefaction
			for (int it = 0; it < 3; ++it) {
				v_r_star.at(it) = (v_r - 2*a_r/(gamma-1)*(1 - pow((p_guess_temp.at(it)/p_r),((gamma-1)/(2*gamma)))));
				rho_r_star.at(it) = (right.rho*pow((p_guess_temp.at(it)/p_r),(1/gamma)));
				E_r_star.at(it) = ((p_guess_temp.at(it)/(gamma-1)) - rho_r_star.at(it)*pow(v_r_star.at(it),2)/2);			
			}
		}
		else { //The right state undergoes a shock
			A = 2/((gamma-1)*right.rho);
			B = (gamma - 1)/(gamma + 1)*p_r;
			for (int it = 0; it < 3; ++it) {
				v_r_star.at(it) = (v_r + (p_guess_temp.at(it)-p_r)*sqrt(A/(p_guess_temp.at(it)+B)));
				rho_r_star.at(it) = (right.rho*((p_r*(gamma-1) + p_guess_temp.at(it)*(gamma+1))/(p_guess_temp.at(it)*(gamma-1) + p_r*(gamma+1))));
				E_r_star.at(it) = ((p_guess_temp.at(it)/(gamma-1)) - rho_r_star.at(it)*pow(v_r_star.at(it),2)/2);
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
			error.at(it) = (v_l_star.at(it) - v_r_star.at(it));
		}

		std::cout << "error = " << error.at(1) << '\n';
/*		for(auto iF:error) {
				std::cout << iF << '\n';
		}
*/

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
			U_star_l.rho = (rho_l_star.at(1)/2); 
			U_star_l.rhou = (U_star_l.rho*v_l_star.at(1));
			U_star_l.E = (p_guess_temp.at(1)/(gamma-1)-(pow(U_star_l.rhou,2)/(U_star_l.rho*2))); 
			U_star_l.pressure = p_guess_temp[1];
			a_star_l = sqrt(gamma*p_guess_temp.at(1)/U_star_l.rho); //sqrt(gamma*P/rho)

			U_star_r.rho = (rho_r_star.at(1)/2); 
			U_star_r.rhou = (U_star_r.rho*v_r_star.at(1));
			U_star_r.E = (p_guess_temp.at(1)/(gamma-1)-(pow(U_star_r.rhou,2)/(U_star_r.rho*2))); 
			U_star_r.pressure = p_guess_temp[1];
			a_star_r = sqrt(gamma*p_guess_temp.at(1)/U_star_r.rho); //sqrt(gamma*P/rho)

 /* Error Checking
			std::cout << '\n' << '\n' << "Yey, we converged!! U Star Left is: " << '\n';
			for(auto ustar:U_star_l) {
				std::cout << ustar << '\n';
			}
			std::cout << "V_l = " << U_star_l.rhou/U_star_l.rho << '\n';
			std::cout << '\n' << "U Star Right is: " << '\n';
			for(auto ustar:U_star_r) {
				std::cout << ustar << '\n';
			}
			std::cout << "V_r = " << U_star_r.rhou/U_star_r.rho << '\n';
End */

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
	
	// Now compute which state is on the x=0 line
	// First, find and order wavespeeds
	//Note - according to paper, need to use modified wavespeeds!(eq.42) - try with these basic here to see if good solutions
	if (type.left == 1) { //The left state undergoes a rarefaction
		wavespeeds.push_back(left.rhou/left.rho - a_l);
		wavespeeds.push_back(U_star_l.rhou/U_star_l.rho - a_star_l);
	}
	else { //The left state undergoes a shock
		wavespeeds.push_back(left.rhou/left.rho - a_l);
	}

	//average of velocities in left and right star state
	wavespeeds.push_back(((U_star_r.rhou/U_star_r.rho) + (U_star_l.rhou/U_star_l.rho))/2); 

	if (type.right == 1) { //The right state undergoes a rarefaction
		wavespeeds.push_back(U_star_r.rhou/U_star_r.rho + a_star_r);
		wavespeeds.push_back(right.rhou/right.rho + a_r);
	}
	else { //The right state undergoes a shock
		wavespeeds.push_back(right.rhou/right.rho + a_r);
	}
	
	return (compute_zerostate(wavespeeds, left, U_star_l, U_star_r, right, gamma, a_l, a_r));	
}




// ------------------------------------------------------------------------------------------

ODstate compute_zerostate(std::vector<double> wavespeeds, ODstate left, ODstate U_star_l, ODstate U_star_r, ODstate right, double gamma, double a_l, double a_r) {

	int left_wavenumber;
	bool computed = 0;
	ODstate in_raref;
	
	for (unsigned int it = 0; it < wavespeeds.size(); ++it) {
		if (wavespeeds.at(it) < 0) {
			continue;
		}
		else if (wavespeeds.at(it) > 0) { // We found the right wave of the 
			computed = 1;
			left_wavenumber = it;
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
			press_l = (left.E + pow(left.rhou,2)/(left.rho*2))*(gamma-1);
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
			in_raref.rhou = in_raref.rho*(2/(gamma+1)*(-a_l + (gamma-1)/2*(right.rhou/right.rho)));

			double press, press_r;
			press_r = (right.E + pow(right.rhou,2)/(right.rho*2))*(gamma-1);
			press = press_r*pow((2/(gamma+1) - (gamma-1)/(a_l*(gamma+1))*(right.rhou/right.rho)),((2*gamma)/(gamma-1)));

			in_raref.E = (press/(gamma-1) - pow(in_raref.rhou,2)/(2*in_raref.rho));
			return(in_raref);	
		}
		else if ((left_wavenumber == 3) && type.left == 2 && type.right == 2) { // left shock and right shock and x=0 is to the right of everything - should never be accessed
			return(right);
		}
		else if ((left_wavenumber == 4) && type.left == 1 && type.right == 2) { // left rarefaction and right shock and x=0 is to the right of everything - should never be accessed
			return(right);
		}
		else if ((left_wavenumber == 4) && type.left == 1 && type.right == 1) { // left rarefaction and right rarefaction and x=0 is inside right rarefaction
			in_raref.rho = right.rho*pow((2/(gamma+1) - (gamma-1)/(a_r*(gamma+1))*(right.rhou/right.rho)),(2/(gamma-1)));
			in_raref.rhou = in_raref.rho*(2/(gamma+1)*(-a_l + (gamma-1)/2*(right.rhou/right.rho)));

			double press, press_r;
			press_r = (right.E + pow(right.rhou,2)/(right.rho*2))*(gamma-1);
			press = press_r*pow((2/(gamma+1) - (gamma-1)/(a_l*(gamma+1))*(right.rhou/right.rho)),((2*gamma)/(gamma-1)));

			in_raref.E = (press/(gamma-1) - pow(in_raref.rhou,2)/(2*in_raref.rho));
			return(in_raref);
		}
		else if ((left_wavenumber == 4) && type.left == 2) { // left shock and x=0 to the right of everything - should never be accessed!
			return(right);
		}
		else if (left_wavenumber == 5) { // x=0 is to the right of everything - should never be accessed!
			return(right);
		}
		else {
			std::cout << '\n' << '\n' << "ERROR: CANNOT FIND WAVE ON X = 0 IN EXACT_RIEMANN_SOLVER.CPP...EXITING" << '\n';
		}
	}
}

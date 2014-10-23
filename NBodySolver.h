#ifndef ODE_H
#define ODE_H

//#include"NBody_functions.h"
#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>

/*
 *      NBodySolver v1.1
 * Uses block timesteps.
 * Solves the equation f'(r_i) = rhs(r_i)
 * for 0 < t < T with minimum timestep dt.
 * f is a 6xN system of first order eqns
 */

using arma::vec;
using arma::mat;
using arma::norm;
using namespace std;

class NBodySolver
{
	private:

		int N;					// #bodies
		double global_t;		// global time
		double T;				// Final time
		double dt;				// constant timestep
		
		mat (*rhs)(mat states, vec masses); // callable object representing the rhs.

		int n_timesteps;
		mat states;
		vec masses;
		vec timesteps;

		
	public:
	
		vector<Body> bodies;
		
		/*
		 * Constructor. Initializes the numerical paramaters
		 */
		
		NBodySolver(int N, mat (*rhs)(mat states, vec masses), double T, double dt);
		
		/*
		* Destructor. Destroy the system
		*/
		
		~NBodySolver();

		/*
		 * Set all the 'bodies' elements w/ masses & initial state
		 */
		
		void setInitialConditions(const char* file);
		
		/*
		 * While global_t is less than final time T; advance the bodies.
		 * This iteratively updates the 'bodies' object.
		 */
		  
		void solve();

		/*
		 * Advance the solutions with Euler's method.
		 */
		
		void advance();
		
		/*
		 * Advance the solutions with RK4 method.
		 */
		
		vec rk4(int i, double dt);
		void recomputeForces();
		void recomputeTimesteps();
		double roundBestTimestep(double dt);
	
		/*
		 * Iterate over the bodies and write them to .csv-file
		 */
		
		void writeBodies(const char * filename);
	
		/*
		 * Writes each body to cout. For developement purposes.
		 */
	
		void writeBodies();
		void writeBodyTrajectory();
	
};

#endif
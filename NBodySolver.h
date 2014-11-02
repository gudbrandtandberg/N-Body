#ifndef ODE_H
#define ODE_H

#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>
#include<set>

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

		int N;
		double global_t;
		double T;
		double G;
		double dtmax;
		double dtmin;
		bool adaptive;
		int n_timesteps;
		mat states;
		vec masses;
		vec timesteps;
		vector<Body> bodies;
		set<int> toStep;
	
	public:
		
		/*
		 * Constructor. Initializes the numerical paramaters
		 */
		
		NBodySolver(int N, double T, double dtmax, bool adaptive);
		
		/*
		* Destructor. Destroy the system
		*/
		
		~NBodySolver();

		/*
		 * Set all the 'bodies' elements w/ masses & initial state
		 */
		
		void setInitialConditions(char* file);
		
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
		 * Advance the solutions with RK4 method. This sets the force attribute on
		 * all bodies in toStep.
		 */
		
		void rk4();
	
		void Verlet();
	
		/*
		 * Just a shell to prepare toStep before calling rk4()
		 */
	
		void recomputeForces();
	
		/*
		 * Checks which bodies needs to recompute their timestep and sets their 
		 * dt and nextEvalTime attributes.
		 */
	
		void recomputeTimesteps();
	
		/*
		 * Returns the forces on the bodies given in the state-matrix states.
		 * Only calculates the forces on the bodies in toStep.
		 */
	
		mat gravity(mat states);
	
		/*
		 * Rounds the number dt to the nearest number in timesteps.
		 */
	
		double roundBestTimestep(double dt);
	
		/*
		 * Iterate over the bodies and write them to .csv-file
		 */
		
		void writeBodies(char * filename);
	
		/*
		 * Writes each body to cout. For developement purposes.
		 */
	
		void writeBodies();
		void writeBodyTrajectory();
	
};

#endif
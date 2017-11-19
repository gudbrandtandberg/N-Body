/*
 * NBodySolver.h
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 * 
 */

#ifndef NBODY_H
#define NBODY_H
#define ARMA_NO_DEBUG

#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>
#include<set>
#include<iomanip>
#include<ctime>
#include<cmath>

#define VERLET 0
#define RK4 1

using namespace arma;
using namespace std;

class NBodySolver
{
	private:

		int N;
		double global_t;
		double T;
		double G;
		double dt;
		double eps;
		mat states;
		mat extrap_states;
		vec masses;
		vec a;
		vec next_state;
		vec r_ij;
		mat rhs;
		mat a_now;
		mat a_next;
		mat v_half;
		mat K1, K2, K3, K4, forces;
		mat positions;
		vector<Body> bodies;
		clock_t start, stop;
		int method;
	
	public:
		
		/*
		 * Constructor. Initializes the numerical paramaters
		 */
		
		NBodySolver(int N, double T, double dt, int method);
		
		/*
		* Destructor. Destroy the system.
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
		 * Advance the solutions with appropriate method.
		 */
		
		void advance();
		
		/*
		 * Advance the solutions with RK4 method. This sets the force attribute on
		 * all bodies in toStep.
		 */
		
		void rk4();
	
		/*
		 * Advance the solutions with Verlet method. This sets the force attribute on
		 * all bodies in toStep.
		 */
	
		void Verlet();
	
		/*
		 * Return the gravity acting on all the bodies.
		 */
	
		mat gravity(mat state);
	
		/*
		 * Iterate over the bodies and write them to .dat-file
		 */
		
		void writeTrajectories();

		/*
		 * Write the total mechanical energy of the system to a file.
		 */
	
		void writeEnergy();
	
};

#endif
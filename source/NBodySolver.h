/*
 * NBodySolver.h
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 * Uses block timesteps. Solves the equation (states)' = gravity(states)
 * for 0 < t < T with maximum timestep dtmax.
 * 
 */

#ifndef ODE_H
#define ODE_H

#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>
#include<set>
#include<iomanip>
#include<ctime>
#include<cmath>
#include<mpi.h>

#define VERLET 0
#define RK4 1

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
		double eps;
		double dtmax;
		double dtmin;
		bool adaptive;
		bool computedFirst;
		int n_timesteps;
		int time_index;
		int mainNode_n;
		double current_dt;
		mat states;
		mat extrap_states;
		vec masses;
		vec a;
		vec r_ij;
		mat rhs;
		mat a_now;
		mat a_next;
		mat v_half;
		mat r_next;
		mat pos;
		vec timesteps;
		mat step_history;
		vector<Body> bodies;
		set<int> toStep;
		clock_t start, stop;
		int method;
	
	public:
		
		/*
		 * Constructor. Initializes the numerical paramaters
		 */
		
		NBodySolver(int N, double T, double dtmax, bool adaptive, int method);
		
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
		 * Just a shell to prepare toStep before calling rk4()
		 */
	
		void recomputeForces();
	
		/*
		 * Checks which bodies needs to recompute their timestep and sets their 
		 * dt and nextEvalTime attributes.
		 */
	
		void computeFirstTimesteps();
	
		void recomputeTimesteps();
	
		/*
		 * Returns the forces on the bodies given in the state-matrix states.
		 * Only calculates the forces on the bodies in toStep.
		 */
	
		mat gravity(mat state);
	
		/*
		 * Rounds the number dt to the nearest number in timesteps.
		 */
	
		double roundBestTimestep(double dt);
	
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
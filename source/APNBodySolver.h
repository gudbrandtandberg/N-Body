/*
 * APNBodySolver.h
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 * Uses block timesteps. Solves the equation (states)' = gravity(states)
 * for 0 < t < T with maximum timestep dtmax.
 *
 */

#ifndef APNBODY_H
#define APNBODY_H

#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>
#include<set>
#include<iomanip>
#include<ctime>
#include<cmath>
//#include<omp.h>

#define VERLET 0
#define RK4 1

using arma::vec;
using arma::mat;
using arma::norm;
using namespace std;

class APNBodySolver
{
private:
	
	int N;
	int cpus;
	
	double global_t;
	double T;
	double G;
	double eps;
	double dtmax, dtmin, dtmed;
	double current_dt;
	
	bool computedFirst;
	
	vec masses;
	vec a;
	vec r_ij;
	vec next_state;
	vec timesteps;
	
	mat extrap_positions;
	mat rhs;
	mat a_now;
	mat a_next;
	mat pos;
	mat step_history;

	vector<Body> bodies;
	
	set<int> min;
	set<int> med;
	set<int> max;
	
	clock_t start, stop;
	
public:
	
	/*
	 * Constructor. Initializes the numerical paramaters
	 */
	
	APNBodySolver(int N, double T, double dtmax, int cpus, double epsilon);
	
	/*
		* Destructor. Destroy the system.
		*/
	
	~APNBodySolver();
	
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
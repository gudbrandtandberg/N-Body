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
#include<omp.h>

#define VERLET 0
#define RK4 1

using arma::vec;
using arma::mat;
using arma::norm;
using namespace std;

class APNBodySolver
{
private:
	
	int N, threads;
	
	double global_t, T, G, eps, eta;
	double dtmax, dtmin, dtmed, current_dt;
	
	bool computedFirst;
	
	vec masses;
	vec a;
	vec next_state;
	vec timesteps;
	vec pos_now;
	vec v_half_now;
	
	mat extrap_positions;
	mat rhs;
	mat a_now;
	mat a_next;
	mat pos;

	vector<Body> bodies;
	
	set<int> min;
	set<int> med;
	set<int> max;
	
	clock_t start, stop;
	
public:
	
	/*
	 * Constructor. Initializes the numerical paramaters
	 */
	
	APNBodySolver(int N, double T, int cpus, double epsilon);
	
	/*
	 * Destructor. Destroy the system.
	 */
	
	~APNBodySolver();
	
	/*
	 * Set all the 'bodies' elements w/ masses & initial state
	 */
	
	void setInitialConditions(char* file);
	
	/*
	 * While global_t is less than final time T: advance the bodies.
	 * This iteratively updates the 'bodies' object.
	 */
	
	void solve();
	
	/*
	 * Computes the acceleration at initial condition and uses this to set a_now and
	 * dt attribute. dtmax is determined by finding the smallest acceleration amongst
	 * the bodies.
	 */
	
	void computeFirstTimesteps();
	
	/*
	 * Resets the dt attribute of each body by using norm(eta/a_now)
	 */
	
	void recomputeTimesteps();
	
	/*
	 * Returns the forces acting on the bodies given in the matrix state.
	 * Only calculates the forces on the bodies with body.dt == current_dt.
	 * Returns zero column for other bodies.
	 */
	
	mat gravity(mat state);
	
	/*
	 * Rounds the number dt to the nearest number in timesteps.
	 */
	
	double roundBestTimestep(double dt);
	
	/*
	 * Iterate over the bodies and write their trajectories to .dat-file
	 */
	
	void writeTrajectories();
	
	
	/*
	 * Do energy/boundedness calculations and write the total mechanical energy, 
	 * the virial energy and the number of bound bodies of the system to a file.
	 */
	
	void writeEnergy();
	
};

#endif
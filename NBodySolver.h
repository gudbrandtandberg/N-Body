#ifndef ODE_H
#define ODE_H

//#include"NBody_functions.h"
#include"Body.h"
#include<fstream>
#include<iostream>
#include<armadillo>

/*
 *      NBodySolver v1.0
 * Uses constant timestep.
 * Solves the equation f'(r_i) = rhs(r_i)
 * for 0 < t < T in n solution steps. 
 * f is a 6xN system of first order eqns
 */

using arma::vec;
using arma::mat;
using arma::norm;
using namespace std;

class NBodySolver;

class NBodySolver
{
	
private:

	int N;					// #bodies
	double global_t;		// global time
	double T;				// Final time
	double dt;				// constant timestep
	
	mat (*rhs)(mat states, vec masses); // callable object representing the rhs.

	vector<Body> bodies;
	vector<Body> toStep;
	mat states;
	vec masses;

	
public:
	
	
	/*
	 * Constructor. Initializes the numerical paramaters
	 */
	
	NBodySolver(int N, mat (*rhs)(mat states, vec masses), double T, double dt);
	
	/*
	* Destructor. Destroy the system
	*/
	
	~NBodySolver();

	/*
	 * Set all the 'bodies' elements w/ masses, initial state, name
	 */
	
	void setInitialConditions(const char* file);
	
	/*
	 * While global_t is less than final time T; advance the bodies.
	 * This iteratively updates the 'bodies' object.
	 */
	  
	void solve();

	/*
	 * Advance the solutions with appropriate method.
	 */
	
	void advanceEuler(double dt);
	
	/*
	 * Iterate over the bodies and write them to .csv-file
	 */
	
	void writeBodies(const char * filename);
	void writeBodies();
	void writeBodyTrajectory();
	
};

#endif
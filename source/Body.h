#ifndef BODY_H
#define BODY_H

#include<iostream>
#include<armadillo>
#include<fstream>

using std::cout;
using std::endl;
using arma::mat;
using arma::vec;
using arma::zeros;

/*
 * Each of the N bodies is represented by a Body object containing information
 * about the mass, times, current state and past states.
 */

class Body
{
public:
	
	int n;
	double dt;
	bool bound;
	
	mat state_history;
	
	vec r, v, v_half, a_now, a_next, force;

	/*
	 * Constructor. Sets initial state and empty trajectory
	 */
	
	Body(vec init_state);
	Body();
	
	/*
	 * Destructor. Destroy the body.
	 */
	
	~Body();
	
	/*
	 * Each time a new state has been calculated by NBodySolver, it is added to the body object.
	 */
	
	void addState(vec state);
	
};

#endif
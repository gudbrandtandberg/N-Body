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
private:
	
	double mass;
	
public:
	
	int n;
	mat state_history;
	vec r;
	vec v;
	vec v_half;
	vec a_now;
	vec a_next;
	double dt;
	vec force;
	bool bound;
	
	/*
	 * Constructor. Sets initial state and empty trajectory
	 */
	
	Body(double mass, vec init_state);
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
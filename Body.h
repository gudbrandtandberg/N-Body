
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
		int n;
		vec t;				 // [t0 t1 t2 ... tn]
		double mass;
		
	public:
		mat trajectory;      // [[x1 y1 z1] [x2 y2 z2] [x3 y3 z3] ... [xt yt zt]]
		vec state;			 // [xt yt zt vxt vyt vzt]
	
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
	
		/*
		 * Print methods for development purposes
		 */
	
		void print();
		void printTrajectory();
	
};

#endif

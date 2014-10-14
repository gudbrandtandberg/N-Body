
#ifndef BODY_H
#define BODY_H

#include<iostream>

#include<armadillo>
#include<fstream>

using std::cout;
using arma::mat;
using arma::vec;
using std::endl;
using arma::zeros;

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
	
		void print();
		void printTrajectory();
		void addState(vec state);
	
};


#endif

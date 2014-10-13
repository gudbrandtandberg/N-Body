
#ifndef BODY_H
#define BODY_H

#include<iostream>

#include<armadillo>
#include<fstream>

using std::cout;
using namespace arma;


class Body
{
	
	private:
		
		vec state;			 // [xt yt zt vxt vyt vzt]
		mat trajectory;      // [[x1 y1 z1] [x2 y2 z2] [x3 y3 z3] ... [xt yt zt]]
		vec t;				 // [t0 t1 t2 ... tn]
		double mass;
		
	public:
		
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
	
};


#endif

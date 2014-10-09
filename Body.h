
#ifndef BODY_H
#define BODY_H

#include<stdio>
#include<ofstream>

using std::cout;
using arma::vec;
using arma::mat;

class Body;

class Body
{
	
private:
	
	vec state;			 // [xt yt zt vxt vyt vzt]
	mat trajectory;      // [[x1 y1 z1] [x2 y2 z2] [x3 y3 z3] ... [xt yt zt]]
	vec t;				 // [t0 t1 t2 ... tn]
	double mass;
	const char * name;
	
public:
	
	/*
	 * Constructor. Sets initial state and empty trajectory
	 */
	
	Body(const char* name, double mass, vec init_state);
	
	/*
	 * Destructor. Destroy the body. 
	 */
	
	~Body();
	
};


#endif

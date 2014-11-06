/*
 * Body.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 */

#include "Body.h"

Body::Body(double mass, vec init_state)
{
	n = 1;
	this->mass = mass;
	this->state = init_state;
	state_history = zeros(6, n);
	state_history.col(0) = state;
	nextEvalTime = 0;
	force = zeros(6);
}

Body::Body()
{
	
}

Body::~Body()
{
	
}

void Body::addState(vec state){

	n = state_history.n_cols;
	state_history.insert_cols(n, state);
	this->state = state;
}

void Body::setNextEvalTime(double time){
	this->nextEvalTime = time;
}

void Body::print()
{
	cout << "======Body=======" << endl;
	cout << "mass: " << mass << endl;
	cout << "state: " <<endl;
	cout << state << endl;
}

void Body::printTrajectory()
{
	for (int i=0; i<state_history.n_cols; i++) {
		cout << state_history.col(i).t() << endl;
	}
}
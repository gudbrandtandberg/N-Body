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
	
	state_history = zeros(6, n);
	state_history.col(0) = init_state;
	a_now = zeros(3);
	a_next = zeros(3);
	v = init_state.rows(3, 5);
	v_half = zeros(3);
	r = init_state.rows(0, 2);
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
	r = state.rows(0, 2);
	v = state.rows(3, 5);
	
}

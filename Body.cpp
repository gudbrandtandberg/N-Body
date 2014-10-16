#include "Body.h"


Body::Body(double mass, vec init_state)
{
	n = 1;
	this->mass = mass;
	this->state = init_state;
	trajectory = mat(3, n);
	trajectory.col(0) = state.rows(0,2);
}

Body::Body()
{
	
}

Body::~Body()
{
	
}

void Body::addState(vec state){
	
	trajectory.insert_cols(n, state.rows(0,2));
	n += 1;
	this->state = state;
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
	for (int i=0; i<n; i++) {
		cout << trajectory.col(i).t() << endl;
	}
}
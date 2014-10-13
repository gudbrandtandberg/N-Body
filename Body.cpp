#include "Body.h"


Body::Body(double mass, vec init_state)
{
	this->mass = mass;
	this->state = init_state;
	
}

Body::Body()
{
		
}

Body::~Body()
{
	
}

void Body::print()
{
	cout << "======Body=======" << endl;
	cout << "mass: " << mass << endl;
	cout << "state: " <<endl;
	cout << state << endl;
}

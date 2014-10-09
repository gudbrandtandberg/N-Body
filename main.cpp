#include"NBody_functions.h"
#include"NBodySolver.h"

/*
* Main program for Gravitational N-body simulation project by Gudbrand Tandberg
* Version 1.0
*
*/


/*
  Questions:
 
 
*/

using arma::vec;
using arma::mat;
using arma::norm;

int main(int argc, char argv[])
{
	
	// some particular test-case:
	int T = 1;
	double dt = 0.1;
	mat (*gravity)(vector<Body>) = gravity;
	int N = 1;
	
	
	NBodySolver solver = NBodySolver(N, gravity, T, dt);
	solver.setInitialConditions("earth_sun.csv");
	solver.solve();
	
}
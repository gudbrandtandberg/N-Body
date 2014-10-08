#include"NBody_functions.h"
#include"ODESolver.h"

/*
* Main program for Gravitational N-body simulation project by Gudbrand Tandberg
* Version 1.0
*
*
*/

/*
 
 Some tough & some not so tough design decisions:
 
 - write solver for both first and second order ODE's? Only need second.
 - write solver for both systems and scalar ODE's? How? Input #eqs and #dims perhaps
 - what kind of structure to store the postisions/velocities in? Keep in mind - mutable.
 - possibility: write a 'Body' class with mass (scalar), positions & velocities (matrices), times (vector)
 - store current state in a 'Body *bodies'. Is this possible in C++?
 - write different methods as subclasses? make ODESolver abstract.
 - what kind of files to write initial conditions? .csv most prob, but to include #eqs & #dims in first line?
 - using function pointer to represent the rhs OK? ptr *gravity
 
 */


int main(int argc, char argv[])
{
	
	// some particular test-case:
	int T = 1;
	double n = 100;
	function_pointer shm = shm;   //implement the rhs. of x'' = x in some function in NBody_functions..
	
	
	ODESolver solver = ODESolver(shm, T, n);
	solver.setMethod(RK4);
	
	solver.setInitialConditions("file.txt");
	solver.solve();
	
	
	double ***POSITIONS = ODESolver.POSITIONS; // further analyze results in matlab
	
	
	
}
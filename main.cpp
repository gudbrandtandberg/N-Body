#include"NBody_functions.h"
#include"NBodySolver.h"
#include<cmath>

/*
* Main program for Gravitational N-body simulation project by Gudbrand Tandberg
* Version 1.0
*
*/

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;

mat gravity(mat states, vec masses)
{
	// use direct summation to calculate force on each body
	
	double G = 1.0;
	int n_bodies = states.n_cols;
	
	mat rhs = zeros(6, n_bodies);
	
	// copy velocities
	rhs.submat(0, 0, 2, n_bodies-1) += states.submat(3, 0, 5, n_bodies-1);
	
	// then determine acceleration
	mat pos = zeros(3, n_bodies);
	mat ax = zeros(3, n_bodies);
	
	pos += states.submat(0, 0, 2, n_bodies-1);
	
	
	vec r_ij = zeros(3);
	vec a = zeros(3);
	
	for (int i=0; i<n_bodies; i++) {  //calc. force on body i
		a = zeros(3);
		
		for (int j=0; j<n_bodies; j++) { //add contribution from all j != i
			if (i != j) {
				
				r_ij = pos.col(j)-pos.col(i);
				a += (G*masses(j)/pow(norm(r_ij), 3))*r_ij;
				
			}
		}
		ax.col(i) = a;
	}
	
	rhs.submat(3, 0, 5, n_bodies-1) = ax;
	
	return rhs;
}


int main(int argc, char **argv)
{
	hei();
	// some particular test-case:
//	int T = 100;
//	double dt = 0.01;
//	mat (*gr)(mat, vec) = gravity;
//	int N = 18;
//
//	
//	NBodySolver solver = NBodySolver(N, gr, T, dt);
//	solver.setInitialConditions("./initial_conditions/solarsystem18.csv");
//	solver.solve();
//	solver.writeBodies("./trajectories/solarsystem18_trajectories.dat");
	
	
}
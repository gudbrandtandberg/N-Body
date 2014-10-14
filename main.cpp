//#include"NBody_functions.h"
#include"NBodySolver.h"
#include<cmath>

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
using arma::zeros;

//G_s = G*d^2/AU^3;           %scaled gravitational constant
//
//m = [1.9891E30 3.302E23 4.8685E24 5.97219E24 6.4185E23 1.8986E27 5.6846E26...
//	 8.681E25 1.0243E26];    %masses of the sun and our 8 planets
//a = zeros(1, planets*3);    %preallocation
//
//%This loop iterates over all the planets adding the acceleration of planet
//%i due to planet j to the overall acceleration on planet i.
//for i = 1:planets
//for j = 1:planets
//if j ~= i
//r_ij = r((j-1)*3+1:j*3) - r((i-1)*3+1:i*3);
//a((i-1)*3+1:i*3) = a((i-1)*3+1:i*3) + (G_s*m(j)/(norm(r_ij)^3))*r_ij;
//end
//end
//end
//end

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
	
	// some particular test-case:
	int T = 1000;
	double dt = 0.05;
	mat (*gr)(mat, vec) = gravity;
	int N = 11;

	
	NBodySolver solver = NBodySolver(N, gr, T, dt);
	solver.setInitialConditions("./initial_conditions/solarsystem11.csv");
	solver.solve();
	solver.writeBodies("./trajectories/solarsystem11_trajectories.dat");
	
	
}
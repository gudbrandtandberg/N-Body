#include"NBody_functions.h"
#include"NBodySolver.h"
#include<cmath>

/*
* Main program for Gravitational N-body simulation project by Gudbrand Tandberg
* Version 1.1
*
*/

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;
using arma::ones;

mat gravity(mat states, vec masses)
{
	// use direct summation to calculate force on each body
	
	double G = 1;
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

//vec gravity(mat states, vec masses, int body)
//{
//	// use direct summation to calculate force on body #i.
//	
//	double G = 1.0;
//	int n_bodies = states.n_cols;
//	
//	vec rhs = zeros(6);
//	
//	// copy velocities
//	rhs.rows(0, 2) += states.col(body).rows(3, 5);
//	
//	// then determine acceleration
//	mat pos = zeros(3, n_bodies);
//	pos += states.submat(0, 0, 2, n_bodies-1);
//	
//	vec r_ij = zeros(3);
//	vec a = zeros(3);
//
//	for (int j=0; j<n_bodies; j++) { //add contribution from all j != i
//		if (j != body) {
//				
//			r_ij = pos.col(j)-pos.col(body);
//			a += (G*masses(j)/pow(norm(r_ij), 3))*r_ij;
//				
//		}
//	}
//	
//	rhs.rows(3, 5) += a;
//	return rhs;
//}



int main(int argc, char **argv)
{
	// some particular test-case:
	int T = 15;
	double dt = 0.01;
	mat (*gr)(mat, vec) = gravity;
	int N = 11;
	
	NBodySolver solver = NBodySolver(N, gr, T, dt);
	solver.setInitialConditions("./initial_conditions/solarsystem11.csv");
	solver.solve();
	solver.writeBodies("./trajectories/solarsystem11_trajectories.dat");
	
	
	
	
//	solver.setInitialConditions("./initial_conditions/3body.csv");
//	solver.solve();
//	solver.writeBodies("./trajectories/3body_trajectories.dat");

//	
//	vec init = ones(6);
//	Body b = Body(1.0, init);
//	
//	b.setNextEvalTime(2);
//	cout << b.nextEvalTime << endl;
//
//	vec sec = init + init;
//	b.addState(sec);
//	b.addState(sec);
//	b.print();
//	b.printTrajectory();
//	
//	cout << b.state_history << endl;
//
//	vector<Body> bodies = vector<Body>();
//	bodies.push_back(b);
//	
//	for (Body bod: bodies){
//		bod.dt = 2.7;
//	}
//	
//	for (int i=0; i<1; i++) {
//		bodies[i].dt = 1;
//	}
//	
//	bodies[0].dt = 1;
//	
//	cout << bodies[0].dt << endl;

}


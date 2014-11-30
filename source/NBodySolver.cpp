/*
 * NBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 */

#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, double T, double dt, int method){
		
	this->T = T;
	this->N = N;
	this->dt = dt;
	this->method = method;
	global_t = 0.0;
	
	if (N >= 100){
		G = pow(M_PI, 2)*pow(20, 3)/(80*N);
		eps = 0.5*10E-1;
	}
	else {
		G = 1.0;
		eps = 0;
	}
	
	bodies = vector<Body>();
	states = zeros(6, N);
	masses = zeros(N);
	positions = zeros(3, N);
	a_now = zeros(3, N);
	a_next = zeros(3, N);
	next_state = zeros(6);
	a = zeros(3);
	r_ij = zeros(3);
	v_half = zeros(6, N);
	rhs = zeros(6, N);
	
}

NBodySolver::~NBodySolver(){
	
	// Destroy the system.
}
	
void NBodySolver::setInitialConditions(char* file){
	
	ifstream infile(file);

	Body b;
	vec init_state = zeros(6);
	double mass, x, y, z, vx, vy, vz;

	int i = 0;
	while (infile >> mass >> x >> y >> z >> vx >> vy >> vz) {
		
		masses(i) = mass;
		init_state << x << y << z << vx << vy << vz;
		
		// create and store a body
		b = Body(init_state);
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
		i++;
	}
	
}

	
void NBodySolver::solve() {

	start = clock();
	
	if (method == VERLET){
		a_now = gravity(states);
		for (int i=0; i<N; i++){
			bodies[i].a_now = a_now.col(i).rows(3, 5);
		}
	}
	
	while (global_t <= T){
		
		printf("Progress: %4.1f %% \r", 100*global_t/T);
		
		switch (method){
			case VERLET:
				Verlet();
				break;
				
			case RK4:
				rk4();
				break;
		}
		
		advance();
		global_t += dt;
	
	}
	stop = clock();
	
	cout << "Solved " << N << "-body system in " << (double)(stop-start)/CLOCKS_PER_SEC << " seconds" << endl;
	cout << "using " << (method == 0 ? "Verlet method" : "RK4 method") << endl;
	cout << "T = " << T << ", N = " << N << " and dt = " << dt << endl;

}


void NBodySolver::advance(){
	
	for (int i=0; i<N; i++) {
		next_state.rows(0, 2) = bodies[i].r + dt*bodies[i].force.rows(0, 2);
		next_state.rows(3, 5) = bodies[i].v + dt*bodies[i].force.rows(3, 5);
		bodies[i].addState(next_state);
		states.col(i) = next_state;
	}
}


void NBodySolver::rk4()
{

	K1 = gravity(states);
	K2 = gravity(states + 0.5*dt*K1);
	K3 = gravity(states + 0.5*dt*K2);
	K4 = gravity(states + dt*K3);
	
	forces = 1/6.0*(K1 + 2*K2 + 2*K3 + K4);
	
	for (int i=0; i<N; i++) {
			bodies[i].force = forces.col(i);
	}
}

void NBodySolver::Verlet()
{
	for (int i=0; i<N; i++) {
		v_half.col(i).rows(0, 2) = bodies[i].v + 0.5*dt*bodies[i].a_now;
	}
	
	a_next = gravity(states + v_half*dt);
	
	for (int i=0; i<N; i++) {
		bodies[i].force.rows(0, 2) = v_half.col(i).rows(0, 2);
		bodies[i].force.rows(3, 5) = (bodies[i].a_now + a_next.col(i).rows(3, 5))/2;
		bodies[i].a_now = a_next.col(i).rows(3, 5);
	}

}

mat NBodySolver::gravity(mat state)
{
	// use direct summation to calculate force on each body in toStep
	
	positions = state.submat(0, 0, 2, N-1);
	
	for (int i=0; i<N; i++) {  //calc. force on body i
			
		a = zeros(3);
			
		for (int j=0; j<N; j++) { //add contribution from all j != i
			if (i != j){
				r_ij = positions.col(j)-positions.col(i);
				a += (G*masses(j)/(pow(norm(r_ij), 3) + norm(r_ij)*pow(eps, 2)))*r_ij;
			}
			
			rhs.col(i).rows(0, 2) = state.col(i).rows(3, 5);
			rhs.col(i).rows(3, 5) = a;
		}
	}
	return rhs;
}

void NBodySolver::writeTrajectories(){
	// iterate over bodies and write to file

	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_trajectories_%.0f_%.2f_0_%d.dat", N, T, dt, method);
	
	ofstream f;
	f.open(filename);
	
	int n = bodies[0].state_history.n_cols;
	mat trajectories = zeros(n, 3*N);
	
	for (int i=0; i<N; i++) {
		
		trajectories.cols(i*3, i*3+2) = bodies[i].state_history.rows(0, 2).t();

	}
	
	f << setprecision(15) << trajectories;
	f.close();
	
	cout << "Wrote trajectories to " << filename << endl;
	
}

void NBodySolver::writeEnergy(){
	
	int n = bodies[0].state_history.n_cols;
	
	vec kinetic = zeros(n);
	vec potential = zeros(n);
	vec total_energy = zeros(n);
	
	// compute kinetic energy
	for (int j=0; j<n; j++){
		for (int i=0; i<N; i++){
			
			kinetic[j] += 0.5*masses(i)*pow(norm(bodies[i].state_history.col(j).rows(3, 5)), 2);
		
		}
	}
	
	
	 // compute potential energy
	double r_ij = 0;
	
	for (int k=0; k<n; k++){      //iterate over time
		for (int i=0; i<N; i++){  //iterate over pairs of bodies
			for (int j=i+1; j<N; j++){
				
				r_ij = norm(bodies[i].state_history.col(k).rows(0, 2) - bodies[j].state_history.col(k).rows(0, 2));
				potential[k] -= masses(i)*masses(j)*G/r_ij;
			}
		}
	}
	
	// write total energy to file
	total_energy = potential + kinetic;
	
	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_energy_%.0f_%.2f_0_%d.dat", N, T, dt, method);
	
	ofstream f;
	f.open(filename);
	f << setprecision(15) << total_energy;  //total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
}
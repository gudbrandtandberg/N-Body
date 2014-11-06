/*
 * NBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 *
 */

#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, double T, double dtmax, bool adaptive){
		
	this->T = T;
	this->N = N;
	this->dtmax = dtmax;
	this->adaptive = adaptive;
	
	global_t = 0.0;
	
	//G = 1.0;
	G = 2*4*pow(M_PI, 2)/(32*N*10);
	
	toStep = set<int>();
	bodies = vector<Body>();
	states = zeros(6, N);
	masses = zeros(N);
	
	if (adaptive){
		n_timesteps = 8;
		timesteps = zeros(n_timesteps);
		timesteps(0) = dtmax;
		
		for (int i=1; i<n_timesteps; i++) {
			timesteps(i) = timesteps(i-1)/2;
		}
		dtmin = timesteps(n_timesteps-1);
	}
	else{
		dtmin = dtmax;
		for (int i=0; i<N; i++){
			toStep.insert(i);
		}
	}
	
}

NBodySolver::~NBodySolver(){
	
	// Destroy the system.
}
	
void NBodySolver::setInitialConditions(char* file){
	
	ifstream infile(file);

	Body b;
	vec init_state = zeros(6, 1);
	double mass, x, y, z, vx, vy, vz;

	int i = 0;
	while (infile >> mass >> x >> y >> z >> vx >> vy >> vz) {
		
		masses(i) = mass;
		init_state << x << y << z << vx << vy << vz;
		
		
		// create and store a body (initially with minimum timestep)
		b = Body(mass, init_state);
		b.dt = this->dtmin;
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
		i++;
	}
	
}

void NBodySolver::setMethod(int method){
	
}
	
void NBodySolver::solve(){

	start = clock();
	while (global_t < T){
		recomputeForces();
		advance();
		if (adaptive){
			recomputeTimesteps();
		}
		global_t += dtmin;

	}
	stop = clock();
	
	cout << "Solved " << N << "-body system in " << (stop-start) << " ticks." << endl;
	if (adaptive){
		cout << "Used block timesteps " << timesteps << endl;
	}
}

/*
 * Does not really recompute the forces, it simply determines whether the forces
 * should be recomputed for the bodies. RK4 calculates forces.
 */

void NBodySolver::recomputeForces()
{
	if (adaptive){
		toStep.clear();
		for (int i=0; i<N; i++) {
		
			if ((bodies[i].nextEvalTime - global_t) < 1.E-8) {
				toStep.insert(i);
			}

		}
	}
	rk4();
}

void NBodySolver::advance(){
	
	vec nextState;
	for (int i=0; i<N; i++) {
		nextState = bodies[i].state + dtmin*bodies[i].force;
		bodies[i].addState(nextState);
		states.col(i) = nextState;
	}
}


void NBodySolver::recomputeTimesteps()
{
	for (int i=0; i<N; i++) {
		if ((bodies[i].nextEvalTime - global_t) < 1.E-8) {
			bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].force.rows(3, 5)));
			bodies[i].setNextEvalTime(global_t + bodies[i].dt);
			
			if (i == 1){
				cout << bodies[i].dt << endl;
			}
		}
	}
}

void NBodySolver::rk4()
{
	
	mat K1 = gravity(states);
	mat K2 = gravity(states + 0.5*dtmin*K1);
	mat K3 = gravity(states + 0.5*dtmin*K2);
	mat K4 = gravity(states + dtmin*K3);
	
	mat forces = 1/6.0*(K1 + 2*K2 + 2*K3 + K4);
	
	for (int i=0; i<N; i++) {
		if (toStep.find(i) != toStep.end()){
			bodies[i].force = forces.col(i);
		} //else do nothing (could be done more neatly)
	}
}

void NBodySolver::Verlet()
{
	mat a0 = gravity(states);
	
	mat a1 = a0.submat(3, 0, 5, N-1);
	mat a2 = gravity(states).submat(3, 0, 5, N-1);
	
	mat forces = zeros(6, N);
	forces.submat(3, 0, 5, N-1) = 0.5*(a1 + a2);
	

	forces.submat(0, 0, 2, N-1) = a0.submat(0, 0, 2, N-1) + 0.5*dtmin*a1;

	for (int i=0; i<N; i++) {
		if (toStep.find(i) != toStep.end()){
			bodies[i].force = forces.col(i);
		} //else do nothing (could be done more neatly)
	}

}

mat NBodySolver::gravity(mat states)
{
	// use direct summation to calculate force on each body in toStep
	double eps = 10E-2;
	mat rhs = zeros(6, N);
	
	mat pos = zeros(3, N);
	pos += states.submat(0, 0, 2, N-1);
	
	vec r_ij = zeros(3);
	vec a = zeros(3);
	
	for (int i=0; i<N; i++) {  //calc. force on body i
		
		if (toStep.find(i) != toStep.end()) {
			
			a = zeros(3);
			
			for (int j=0; j<N; j++) { //add contribution from all j != i
				if (i != j) {
					
					r_ij = pos.col(j)-pos.col(i);
					a += (G*masses(j)/(pow(norm(r_ij), 3) + norm(r_ij)*pow(eps, 2)))*r_ij;
					
				}
			}
			rhs.col(i).rows(0, 2) = states.col(i).rows(3, 5);
			rhs.col(i).rows(3, 5) = a;
		}
		else{
			rhs.col(i) = bodies[i].force;
		}
	}
	
	return rhs;
}


double NBodySolver::roundBestTimestep(double dt)
{
	//round dt to block timestep
	double cmp;
	double smallest = abs(dt-timesteps(0));
	int smallestInd = 0;
	
	for (int i=1; i<n_timesteps; i++) {
		
		cmp = abs(timesteps(i)-dt);
		
		if (cmp < smallest){
			smallestInd = i;
			smallest = cmp;
		}
		
	}
	
	return timesteps(smallestInd);
	
}

void NBodySolver::writeTrajectories(){
	// iterate over bodies and write to file

	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_trajectories_%.0f_%.1f_%d.dat", N, T, dtmax, adaptive);
	
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
	
	for (int j=0; j<n; j++){			//iterates over time
		for (int i=0; i<N; i++){		//iterates over bodies
			
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
	sprintf(filename, "./output/%d_body_energy_%.0f_%.1f_%d.dat", N, T, dtmax, adaptive);
	
	ofstream f;
	f.open(filename);
	f << setprecision(15) << total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
}





/*
 * NBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 */

#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, double T, double dtmax, bool adaptive, int method){
		
	this->T = T;
	this->N = N;
	this->dtmax = dtmax;
	this->adaptive = adaptive;
	this->method = method;
	global_t = 0.0;
	
	if (N >= 100){
		G = pow(M_PI, 2)/(2*N*10);
		eps = 10E-2;
	}
	else {
		G = 1.0;
		eps = 0;
	}
	
	toStep = set<int>();
	bodies = vector<Body>();
	states = zeros(6, N);
	masses = zeros(N);
	pos = zeros(3, N);
	a_now = zeros(3, N);
	a_next = zeros(3, N);
	r_next = zeros(3, N);
	computedFirst = false;
	a = zeros(3);
	r_ij = zeros(3);
	v_half = zeros(3, N);
	rhs = zeros(6, N);
	time_index = 0;
	current_dt = 0;
	mainNode_n = 0;
	
	if (adaptive){
		n_timesteps = 3;
		timesteps = zeros(n_timesteps);
		timesteps(0) = dtmax;
		
		for (int i=1; i<n_timesteps; i++) {
			timesteps(i) = timesteps(i-1)/2;
		}
		dtmin = timesteps(n_timesteps-1);
		
		step_history = zeros(N, ceil(T/dtmin)+1);

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
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
		i++;
	}
	
}

	
void NBodySolver::solve() {

	start = clock();
	computeFirstTimesteps();
	while (global_t <= T){
		
		cout << bodies[0].nextEvalTime << endl;
		
		recomputeForces();
		advance();
		if (adaptive && time_index % (int) pow(2, n_timesteps-1) == 0 && time_index != 0){
			recomputeTimesteps();
			mainNode_n++;
		}
		global_t += dtmin;
		time_index++;
		

	}
	stop = clock();
	
	cout << "Solved " << N << "-body system in " << (stop-start) << " ticks" << endl;
	cout << "using " << (method == 0 ? "Verlet method" : "RK4 method") << endl;
	cout << "T = " << T << ", N = " << N << " and dt = " << dtmax << endl;
	if (adaptive){
		cout << "Used block timesteps " << endl << timesteps << endl;
	}

}

void NBodySolver::computeFirstTimesteps()
{

	for (int i=0; i<N; i++){ //initially everyone steps
		toStep.insert(i);
	}
	
	mat forces = gravity(states);
	computedFirst = true;
	a_now = forces.submat(3, 0, 5, N-1);
	
	//cout << a_now << endl;
	
	for (int i=0; i<N; i++) {
			bodies[i].force = forces.col(i);
			bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].force.rows(3, 5)));
			bodies[i].setNextEvalTime(global_t + bodies[i].dt);
		

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
	if (method == VERLET){
		Verlet();
	}
	else{
		rk4();
	}
	
}

void NBodySolver::advance(){
	
	vec nextState;
	for (int i=0; i<N; i++) {
		nextState = bodies[i].state + dtmin*bodies[i].force;
		bodies[i].addState(nextState);
		states.col(i) = nextState;
		step_history(i, time_index) = bodies[i].dt;
	}
}


void NBodySolver::recomputeTimesteps()
{
	for (int i=0; i<N; i++) {
		bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].force.rows(3, 5)));
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
		}
	}
}

void NBodySolver::Verlet()
{

}

mat NBodySolver::gravity(mat state)
{
	// use direct summation to calculate force on each body in toStep
	
	rhs = zeros(6, N);
	pos = state.submat(0, 0, 2, N-1);
	
	
	for (int i=0; i<N; i++) {  //calc. force on body i
		
		if ((toStep.find(i) != toStep.end() && bodies[i].dt == current_dt) || !computedFirst) {
			
			a = zeros(3);
			
			for (int j=0; j<N; j++) { //add contribution from all j != i
				
				if (i != j){
					
					r_ij = pos.col(j)-pos.col(i);
					a += (G*masses(j)/(pow(norm(r_ij), 3) + norm(r_ij)*pow(eps, 2)))*r_ij;
				
				}
			
			rhs.col(i).rows(0, 2) = state.col(i).rows(3, 5);
			rhs.col(i).rows(3, 5) = a;
			//cout << "computed body " << i << endl;
			}
		}
		else {
			rhs.col(i) = zeros(6);
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
	sprintf(filename, "./output/%d_body_trajectories_%.0f_%.1f_%d_%d.dat", N, T, dtmax, adaptive, method);
	
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
		//cout << " Kinetic: "<< kinetic[j] << endl;
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
		//cout << " Potential: "<< potential[k] << endl;
	}
	
	
	// write total energy to file
	total_energy = potential + kinetic;
	
	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_energy_%.0f_%.1f_%d_%d.dat", N, T, dtmax, adaptive, method);
	
	ofstream f;
	f.open(filename);
	f << setprecision(15) << step_history.t();  //total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
}





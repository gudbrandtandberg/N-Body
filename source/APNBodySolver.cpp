/*
 * APNBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 */

#include"APNBodySolver.h"

APNBodySolver::APNBodySolver(int N, double T, double dtmax){
	
	this->T = T;
	this->N = N;
	this->dtmax = dtmax;
	
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
	min = set<int>();
	med = set<int>();
	max = set<int>();
	
	bodies = vector<Body>();
	masses = zeros(N);
	pos = zeros(3, N);
	a_now = zeros(3, N);
	a_next = zeros(3, N);
	computedFirst = false;
	a = zeros(3);
	r_ij = zeros(3);
	next_state = zeros(6);
	rhs = zeros(3, N);
	time_index = 0;
	current_dt = 0;
	extrap_positions = zeros(3, N);
	
	
	n_timesteps = 3;
	timesteps = zeros(n_timesteps);
	timesteps(0) = dtmax;
	
	for (int i=1; i<n_timesteps; i++) {
		timesteps(i) = timesteps(i-1)/2;
	}
	
	dtmin = timesteps(n_timesteps-1);
	dtmed = timesteps(1);
	
	step_history = zeros(N, ceil(T/dtmin)+1);
	
	
}

APNBodySolver::~APNBodySolver(){
	
	// Destroy the system.
}

void APNBodySolver::setInitialConditions(char* file){
	
	ifstream infile(file);
	
	Body b;
	vec init_state = zeros(6);
	double mass, x, y, z, vx, vy, vz;
	
	int i = 0;
	while (infile >> mass >> x >> y >> z >> vx >> vy >> vz) {
		
		masses(i) = mass;
		init_state << x << y << z << vx << vy << vz;
		// create and store a body (initially with minimum timestep)
		b = Body(mass, init_state);
		bodies.push_back(b);
		// also store the initial states in states matrix
		//states.col(i) = init_state;
		i++;
	}
	
}


void APNBodySolver::solve() {
	
	
	start = clock();
	computeFirstTimesteps();  //computes first timesteps and sets a_now, v, r_now for all bodies
	
	int step = 0;
	
	while (global_t < T){
		step = 0; //set to zero before each main step.
		
		for (int i1 = 0; i1<2; i1++) {
			for (int i2=0; i2<2; i2++) {
				step++;
				current_dt = dtmin;
				//update planets with dt = dtmin
				
				for (int i=0; i<N; i++) {
					if (min.find(i) != min.end()) {
						
						bodies[i].v_half = bodies[i].v + bodies[i].a_now*dtmin/2;
						extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*dtmin;
						
					}
					
					else if (med.find(i) != med.end()){
						bodies[i].v_half = bodies[i].v + bodies[i].a_now*dtmed/2;
						extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*(step%2==1?1:2)*dtmin;
					}
					else {
						bodies[i].v_half = bodies[i].v + bodies[i].v*dtmax/2;
						extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*step*dtmin;
					}
					
					
				}
				
				a_next = gravity(extrap_positions);
				
				//walk these planets one step
				for (int i=0; i<N; i++) {
					if (min.find(i) != min.end()) {
						next_state.rows(0, 2) = bodies[i].r + dtmin*bodies[i].v_half;
						next_state.rows(3, 5) = bodies[i].v + dtmin*(bodies[i].a_now + a_next.col(i))/2;
						bodies[i].addState(next_state);
						bodies[i].a_now = a_next.col(i);
					}
				}
				
				global_t += dtmin;
			}
			
			//update planets with dt = dtmed
			current_dt = dtmed;
			
			for (int i=0; i<N; i++) {
				if (med.find(i) != med.end()) {
					
					extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*dtmed;
					
				}
				
				else if (min.find(i) != min.end()){
					extrap_positions.col(i) = bodies[i].r;
				}
				else {
					extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*step*dtmin;
				}
				
				
			}
			
			a_next = gravity(extrap_positions);
			
			//walk these planets one step (twice!)
			for (int i=0; i<N; i++) {
				if (med.find(i) != med.end()) {
					vec pos_now = bodies[i].r;
					vec v_half_now = bodies[i].v;
					for (int j=1; j<3; j++){
						next_state.rows(0, 2) = pos_now + j*dtmin*bodies[i].v_half;
						next_state.rows(3, 5) = v_half_now + dtmed*(bodies[i].a_now + a_next.col(i))/2;
						bodies[i].addState(next_state);
					}
					bodies[i].a_now = a_next.col(i);
				}
			}
		}
		
		//update planets with dt = dtmax
		current_dt = dtmax;
		
		for (int i=0; i<N; i++) {
			if (max.find(i) != max.end()) {
				extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*dtmax;
			}
			
			else {
				extrap_positions.col(i) = bodies[i].r;
			}
		}
		
		a_next = gravity(extrap_positions);
		
		//walk these planets one step (four times!)
		for (int i=0; i<N; i++) {
			if (max.find(i) != max.end()) {
				vec pos_now = bodies[i].r;
				vec v_half_now = bodies[i].v;
				for (int j=1; j<5; j++){
					next_state.rows(0, 2) = pos_now + j*dtmin*v_half_now;
					next_state.rows(3, 5) = v_half_now + dtmax*(bodies[i].a_now + a_next.col(i))/2;
					bodies[i].addState(next_state);
				}
				bodies[i].a_now = a_next.col(i);
			}
		}
		
		recomputeTimesteps(); //at every main node
	}
	stop = clock();
	cout << "Solved " << N << "-body system in " << (stop-start) << " ticks" << endl;
	cout << "using " << "Verlet method" << endl;
	cout << "T = " << T << ", N = " << N << " and dt = " << dtmax << endl;
	cout << "Used block timesteps " << endl << timesteps << endl;
	
	
}

void APNBodySolver::computeFirstTimesteps()
{
	
	mat states = zeros(3, N);
	for (int i=0; i<N; i++) {
		states.col(i) = bodies[i].r;
	}
	mat accelerations = gravity(states);
	computedFirst = true;
	for (int i=0; i<N; i++) {
		bodies[i].a_now = accelerations.col(i);
		bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].a_now));
		cout << 0.001/norm(bodies[i].a_now) << endl;
		if (bodies[i].dt == dtmin){
			min.insert(i);
		}
		else if (bodies[i].dt== dtmed){
			med.insert(i);
		}
		else{
			max.insert(i);
		}
	}
}


void APNBodySolver::recomputeTimesteps()
{
	min.clear(); med.clear(); max.clear();
	cout << "Recomputing " << endl;
	for (int i=0; i<N; i++) {
		bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].a_now));
		cout << 0.001/norm(bodies[i].a_now) << endl;
		if (bodies[i].dt == dtmin){
			cout << i << " min" << endl;
			min.insert(i);
		}
		else if (bodies[i].dt== dtmed){
			cout << i << " med" << endl;
			med.insert(i);
		}
		else{
			cout << i << " max" << endl;
			max.insert(i);
		}
	}
}


mat APNBodySolver::gravity(mat positions)
{
	// use direct summation to calculate force on each body with correct timestep
	
	for (int i=0; i<N; i++) {  //calc. force on body i
		
		if (bodies[i].dt == current_dt || !computedFirst) {
			
			a = zeros(3);
			
			for (int j=0; j<N; j++) { //add contribution from all j != i
				
				if (i != j){
					
					r_ij = positions.col(j)-positions.col(i);
					a += (G*masses(j)/(pow(norm(r_ij), 3) + norm(r_ij)*pow(eps, 2)))*r_ij;
					
				}
				
				rhs.col(i) = a;
			}
		}
	}
	
	return rhs;
}


double APNBodySolver::roundBestTimestep(double dt)
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

void APNBodySolver::writeTrajectories(){
	// iterate over bodies and write to file
	
	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_trajectories_%.0f_%.1f_1_0.dat", N, T, dtmax);
	
	ofstream f;
	f.open(filename);
	
	int n = bodies[1].state_history.n_cols;
	cout << n << endl;
	mat trajectories = zeros(n, 3*N);
	
	for (int j=0; j<n; j++){
		for (int i=0; i<N; i++) {
			
			trajectories.cols(i*3, i*3+2) = bodies[i].state_history.rows(0, 2).t();
			
		}
	}
	
	f << setprecision(15) << trajectories;
	f.close();
	
	cout << "Wrote trajectories to " << filename << endl;
	
}

void APNBodySolver::writeEnergy(){
	
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
	sprintf(filename, "./output/%d_body_energy_%.0f_%.1f_1_0.dat", N, T, dtmax);
	
	ofstream f;
	f.open(filename);
	f << setprecision(15) << total_energy;  //total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
}

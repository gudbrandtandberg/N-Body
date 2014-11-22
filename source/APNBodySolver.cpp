/*
 * APNBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 */

#include"APNBodySolver.h"

APNBodySolver::APNBodySolver(int N, double T, double dtmax, int cpus, double epsilon){
	
	this->N = N;
	this->T = T;
	this->dtmax = dtmax;
	this->cpus = cpus;
	
	current_dt = 0;
	computedFirst = false;
	
	global_t = 0.0;
	
	if (N >= 100){
		G = (pow(M_PI, 2)*pow(20, 3))/(80*N);
		eps = epsilon;
	}
	else {
		G = 1.0;
		eps = 0;
	}
	
	min = set<int>();
	med = set<int>();
	max = set<int>();
	
	bodies = vector<Body>();
	masses = zeros(N);
	
	rhs = zeros(3, N);
	a_now = zeros(3, N);
	a_next = zeros(3, N);
	
	a = zeros(3);
	r_ij = zeros(3);
	next_state = zeros(6);
	extrap_positions = zeros(3, N);
	
	dtmed = dtmax/2;
	dtmin = dtmax/4;
	
	timesteps << dtmin << dtmed << dtmax;
	
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
		b = Body(mass, init_state);
		bodies.push_back(b);
		i++;
	}

	infile.close();
	
}


void APNBodySolver::solve() {
	
	
	start = clock();
	
	computeFirstTimesteps();  //computes first timesteps and sets a_now for all bodies
	
	int step = 0;
	
	while (global_t < T){
		step = 0; //keeps track of where we are within a main-step
		
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
					
					else if (med.find(i) != med.end()){
						bodies[i].v_half = bodies[i].v + bodies[i].a_now*dtmed/2;
						extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*(step%2==1?1:2)*dtmin;
					}
					else {
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
				printf("Progress: %.2f %% \r", 100.0*global_t/T);
				global_t += dtmin;
			}
			
			//update planets with dt = dtmed
			current_dt = dtmed;
			
			for (int i=0; i<N; i++) {
				if (med.find(i) != med.end()) {
					
					extrap_positions.col(i) = bodies[i].r + bodies[i].v_half*dtmed;
					
				}
				
				else if (min.find(i) != min.end()){
					extrap_positions.col(i) = bodies[i].r;
				}
				else {
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
	cout << endl << "Solved " << N << "-body system in " << (stop-start) << " ticks" << endl;
	cout << "using " << "Verlet method" << endl;
	cout << "T = " << T << ", N = " << N << " and dt = " << dtmax << endl;
	cout << "Used block timesteps " << endl << timesteps << endl;
	
	
}

void APNBodySolver::computeFirstTimesteps()
{
	
	mat positions = zeros(3, N);
	
	for (int i=0; i<N; i++) {
		positions.col(i) = bodies[i].r;
	}
	
	a_now = gravity(positions);
	
	computedFirst = true;
	for (int i=0; i<N; i++) {
		bodies[i].a_now = a_now.col(i);
		bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].a_now));

		if (bodies[i].dt == dtmin){
			min.insert(i);
		}
		if (bodies[i].dt == dtmed){
			med.insert(i);
		}
		else {
			max.insert(i);
		}
	}
}


void APNBodySolver::recomputeTimesteps()
{
	min.clear(); med.clear(); max.clear();

	for (int i=0; i<N; i++) {
		bodies[i].dt = roundBestTimestep(0.001/norm(bodies[i].a_now));
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


mat APNBodySolver::gravity(mat positions)
{
	// use direct summation to calculate force on each body with correct timestep
	
	
	for (int i=0; i<N; i++) {  //calc. force on body i
		
		if (bodies[i].dt == current_dt || !computedFirst) {
			
			a = zeros(3);
			
			//#define NUMBER_OF_THREADS 4
			//#pragma omp parallel for private(r_ij) num_threads(NUMBER_OF_THREADS)
			
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
	
	for (int i=1; i<3; i++) {
		
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
	sprintf(filename, "./output/%d_body_trajectories_%.0f_%.3f_1_0_%d_%.2f.dat", N, T, dtmax, cpus, eps);
	
	ofstream f;
	f.open(filename);
	
	int n = bodies[0].state_history.n_cols;

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
	
	double KE = 0;
	double PE = 0;
	double r_ij = 0;
	
	double total_kinetic = 0;
	double total_potential = 0;
	
	vec total_energy = zeros(n);
	vec virial_energy = zeros(n);
	
	vec body_potential = zeros(N);
	vec body_kinetic = zeros(N);
	
	int numbound_now = 0;
	vec numbound = zeros(n);
	
	for (int k=0; k<n; k++){   //iterates over time
		
		// compute total and individual potential energy at time index k
		for (int i=0; i<N; i++){
			for (int j=i+1; j<N; j++){
				
				r_ij = norm(bodies[i].state_history.col(k).rows(0, 2) - bodies[j].state_history.col(k).rows(0, 2));
				
				PE = masses(i)*masses(j)*G/r_ij;
				
				total_potential -= PE;
				
				body_potential[i] -= PE/2;
				body_potential[j] -= PE/2;
			}
		}
		
		//compute kinetic energy of each body and test if body is bound
		numbound_now = N;
		for (int i=0; i<N; i++){
			
			KE = 0.5*masses(i)*pow(norm(bodies[i].state_history.col(k).rows(3, 5)), 2);
			body_kinetic[i] = KE;
			
			if (body_kinetic[i] + body_potential[i] > 0) {
				bodies[i].bound = false;
				numbound_now--;
			}
			
		}
		numbound[k] = numbound_now;
		
		for (int i=0; i<N; i++) {
			if (bodies[i].bound) {
				virial_energy[k] += 2*body_kinetic[i];
				virial_energy[k] += body_potential[i];
			}
		}
		
		virial_energy[k] /= numbound[k];  //ensemble average
		
		total_energy[k] = total_potential + total_kinetic;
		body_potential = zeros(N);
		body_potential = zeros(N);
	}
	
	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_energy_%.0f_%.3f_1_0_%d_%.2f.dat", N, T, dtmax, cpus, eps);
	
	ofstream f;
	f.open(filename);
	f << setprecision(15) << total_energy;  //total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
	sprintf(filename, "./output/%d_body_virialenergy_%.0f_%.3f_1_0_%d_%.2f.dat", N, T, dtmax, cpus, eps);
	f.open(filename);
	f << setprecision(15) << virial_energy;  //total_energy;
	f.close();
	cout << "Wrote virial energy to " << filename << endl;
	
	sprintf(filename, "./output/%d_body_bound_%.0f_%.3f_1_0_%d_%.2f.dat", N, T, dtmax, cpus, eps);
	f.open(filename);
	f << setprecision(15) << numbound;  //total_energy;
	f.close();
	cout << "Wrote number of bound bodies to " << filename << endl;
	
}

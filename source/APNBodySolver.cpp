/*
 * APNBodySolver.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 fall 2014
 */

#include"APNBodySolver.h"

APNBodySolver::APNBodySolver(int N, double T, int threads, double epsilon){
	
	this->N = N;
	this->T = T;
	this->threads = threads;
	
	current_dt = 0;
	computedFirst = false;
	
	global_t = 0.0;
	
	if (N >= 100){
		G = (100*pow(M_PI, 2))/N;
		eps = epsilon;
		eta = 0.1;
	}
	else {
		G = 1.0;
		eps = 0;
		eta = 0.001;
	}
	
	min = set<int>();
	med = set<int>();
	max = set<int>();
	
	bodies = vector<Body>();
	masses = zeros(N);
	
	rhs = zeros(3, N);
	a_now = zeros(3, N);
	a_next = zeros(3, N);
	extrap_positions = zeros(3, N);

	v_half_now = zeros(3);
	pos_now = zeros(3);
	next_state = zeros(6);
	
	dtmax = 0;
	dtmed = 0;
	dtmin = 0;	
	
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
		b = Body(init_state);
		bodies.push_back(b);
		i++;
	}

	infile.close();
	
}


void APNBodySolver::solve() {
	
	
	start = clock();
	computeFirstTimesteps();
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
					pos_now = bodies[i].r;
					v_half_now = bodies[i].v;
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
	
	cout << endl << "Solved " << N << "-body system in " << ((double)stop-start)/(CLOCKS_PER_SEC*threads) << " seconds" << endl;
	cout << "using " << "Verlet method" << endl;
	cout << "using " << threads << " threads" << endl;
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
	
	double smallest_a = norm(a_now.col(1));
	double cmp = 0;
	
	for (int i=2; i<N; i++) {
		cmp = norm(a_now.col(i));
		if (cmp < smallest_a) {
			smallest_a = cmp;
		}
	}

	dtmax = eta/smallest_a;
	dtmed = dtmax/2;
	dtmin = dtmax/4;
	
	timesteps << dtmax << dtmed << dtmin;
	computedFirst = true;
	
	for (int i=0; i<N; i++) {
		bodies[i].a_now = a_now.col(i);
		bodies[i].dt = roundBestTimestep(eta/norm(bodies[i].a_now));

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
		bodies[i].dt = roundBestTimestep(eta/norm(bodies[i].a_now));
		if (bodies[i].dt == dtmin){
			min.insert(i);
		}
		else if (bodies[i].dt == dtmed){
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
			
			double ax = 0;
			double ay = 0;
			double az = 0;
			
#pragma omp parallel for reduction(+:ax, ay, az) num_threads(threads)
			
			for (int j=0; j<N; j++) { //add contribution from all j != i
				
				if (i == j) continue;
				
				double dx = positions(0,j) - positions(0,i);
				double dy = positions(1,j) - positions(1,i);
				double dz = positions(2,j) - positions(2,i);
				
				double dr2 = dx*dx + dy*dy + dz*dz;
				double dr = sqrt(dr2);
				
				double force = G*masses(j)/(dr2*dr + dr*eps*eps);
				
				ax += force*dx;
				ay += force*dy;
				az += force*dz;
			}
			
			rhs(0, i) = ax;
			rhs(1, i) = ay;
			rhs(2, i) = az;
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
	sprintf(filename, "./output/%d_body_trajectories_%.0f_1_0_%d_%.2f.dat", N, T, threads, eps);
	ofstream f;
	f.open(filename);

	int n = bodies[1].state_history.n_cols-4;
	mat trajectories = zeros(n, 3*N);
	
	for (int i=0; i<N; i++) {
		//Costly transposition operation..
		trajectories.cols(i*3, i*3+2) = bodies[i].state_history.rows(0, 2).cols(0, n-1).t();
	}
	
	f << setprecision(15) << trajectories;
	
	f.close();
	
	cout << "Wrote trajectories to " << filename << endl;
	
}

void APNBodySolver::writeEnergy(){
	
	int n = bodies[1].state_history.n_cols-4;
	
	double KE = 0;
	double PE = 0;
	double dr = 0;
	
	double total_kinetic = 0;
	double total_potential = 0;
	double ejected_energy = 0;
	double unbound_final_energy = 0;
	double v = 0;
	
	vec total_energy = zeros(n);
	vec virial_energy = zeros(n);
	
	vec body_potential = zeros(N);
	vec body_kinetic = zeros(N);
	
	int numbound_now = 0;
	vec numbound = zeros(n);
	
	for (int k=0; k<n; k++){   //iterates over time
		
		// compute total and individual potential energy at time k
		for (int i=0; i<N; i++){
			for (int j=i+1; j<N; j++){
				
				dr = norm(bodies[i].state_history.col(k).rows(0, 2) - bodies[j].state_history.col(k).rows(0, 2));
				
				PE = masses(i)*masses(j)*G/dr;
				
				total_potential -= PE;
				
				body_potential[i] -= PE/2;
				body_potential[j] -= PE/2;
			}
		}
		
		//compute kinetic energy of each body and test if body is bound
		numbound_now = N;
		for (int i=0; i<N; i++){
			v = norm(bodies[i].state_history.col(k).rows(3, 5));
			KE = 0.5*masses(i)*v*v;
			body_kinetic[i] = KE;
			total_kinetic += KE;
			
			if (body_kinetic[i] + body_potential[i] > 0) {
				bodies[i].bound = false;
				numbound_now--;
			}
			else {
				bodies[i].bound = true;
			}
			
		}
		numbound[k] = numbound_now;
		
		for (int i=0; i<N; i++) {
			if (bodies[i].bound) {
				virial_energy[k] += 2*body_kinetic[i] + body_potential[i];
			}
		}
		
		virial_energy[k] /= numbound[k];  //ensemble average
		
		// Check how much energy has been ejected at final time.
		
		if (k == n-1){
			for (int i=0; i<N; i++) {
				if (!bodies[i].bound) {
					unbound_final_energy += body_potential[i];
					unbound_final_energy += body_kinetic[i];
				}
			}
			ejected_energy = 100*unbound_final_energy/(total_kinetic + total_potential);
		}
		
		total_energy[k] = total_potential + total_kinetic;
		body_potential = zeros(N);
		body_potential = zeros(N);
		total_kinetic = 0;
		total_potential = 0;
	}
	
	char *filename = new char[50];
	sprintf(filename, "./output/%d_body_energy_%.0f_1_0_%d_%.2f.dat", N, T, threads, eps);
	ofstream f;
	
	f.open(filename);
	f << setprecision(15) << total_energy;
	f.close();
	cout << "Wrote energy to " << filename << endl;
	
	sprintf(filename, "./output/%d_body_virialenergy_%.0f_1_0_%d_%.2f.dat", N, T, threads, eps);
	f.open(filename);
	f << setprecision(15) << virial_energy;  //total_energy;
	f.close();
	cout << "Wrote virial energy to " << filename << endl;
	
	sprintf(filename, "./output/%d_body_bound_%.0f_1_0_%d_%.2f.dat", N, T, threads, eps);
	f.open(filename);
	f << setprecision(15) << numbound;
	f.close();
	cout << "Wrote number of bound bodies to " << filename << endl;
	
	cout << "Roughly " << ejected_energy << "% of the energy has been ejected at final state" << endl;
}

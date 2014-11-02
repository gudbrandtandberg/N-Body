#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, double T, double dtmax, bool adaptive){
		
	this->T = T;
	this->N = N;
	this->dtmax = dtmax;
	this->adaptive = adaptive;
	
	global_t = 0.0;
	G = 1.0;
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
	
	ifstream f;
	f.open(file);
	char *line = new char[100];
	Body b;
	vec init_state = zeros(6, 1);
	double mass;
	char * word;
	
	for (int i=0; i<N; i++){
		
		// read line i -> corresponds to body i
		f.getline(line, 100);
		
		// read line word for word. first mass
		word = strtok(line, ",");
		mass = atof(word);
		masses(i) = mass;
		
		// then initial state
		for (int j=0; j<6; j++) {
			word = strtok(NULL, ",");
			init_state(j) = atof(word);
		}
		
		// create and store a body (initially with minimum timestep)
		b = Body(mass, init_state);
		b.dt = this->dtmin;
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
	}
	
}
	
void NBodySolver::solve(){

	while (global_t < T){
		recomputeForces();
		advance();
		if (adaptive){
			recomputeTimesteps();
		}
		global_t += dtmin;

	}
	
	cout << timesteps << endl;
	
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
					a += (G*masses(j)/pow(norm(r_ij), 3))*r_ij;
					
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

void NBodySolver::writeBodies(char * filename){
	//iterate over bodies and write to file
	
	ofstream f;
	f.open(filename);
	char *buffer = new char[100];
	
	for (int body=0; body<N; body++) {
			
		f << bodies[body].state_history;

	}
	
	f.close();
	
}

void NBodySolver::writeBodies(){
	for (Body b: bodies){
		b.print();
		//cout << b.state_history;
	}
}

void NBodySolver::writeBodyTrajectory()
{
	bodies[1].printTrajectory();
}




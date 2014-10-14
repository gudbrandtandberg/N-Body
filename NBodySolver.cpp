#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, mat (*rhs)(mat states, vec masses), double T, double dt){
		
	this->T = T;
	this->N = N;
	this->rhs = rhs;
	this->dt = dt;
	
	global_t = 0.0;
	
	bodies = vector<Body>();
	toStep = vector<Body>();
	states = zeros(6, N);
	this->masses = zeros(N);
		
}

NBodySolver::~NBodySolver(){
	
	// Destroy the system.
}
	
void NBodySolver::setInitialConditions(const char* file){
	
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
		
		// create and store a body
		b = Body(mass, init_state);
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
	}
	
}
	
void NBodySolver::solve(){
		
	while (global_t < T){
	// determine timesteps
	// quantize timesteps
	// determine which bodies to advance, add them to toStep
		advanceEuler(dt);
		global_t += dt;
	}
	
}
	
void NBodySolver::advanceEuler(double dt){
	
	// create matrix of states = [[x y z vx vy vz] ... [x y z vx vy vz]]
	// take a Euler-step
	// update the bodies in toStep	
	
	//read states form bodies
	
	for (int i=0; i<N; i++){
		states.col(i) = bodies[i].state;
	}
	
	//cout << rhs(states, masses);
	
	states = states + dt*rhs(states, masses);
	
	
	for (int i=0; i<N; i++){
		 bodies[i].addState(states.col(i));
	}
	
}

/*
 * Writes bodies to csv file. Requires constant timestep.
 */
void NBodySolver::writeBodies(const char * filename){
	//iterate over bodies and write to file
	
	ofstream f;
	f.open(filename);
	char *buffer = new char[100];
	
	for (int body=0; body<N; body++) {
			
		f << bodies[body].trajectory;
		
	}
	
	f.close();
	
}

/*
 * Writes each body to cout. For developement purposes.
 */
void NBodySolver::writeBodies(){
	for (Body b: bodies){
		b.print();
	}
}

void NBodySolver::writeBodyTrajectory()
{
	bodies[1].printTrajectory();
}




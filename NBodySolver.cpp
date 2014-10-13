#include"NBodySolver.h"
#include<string>
	
NBodySolver::NBodySolver(int N, mat (*rhs)(vector<Body>), double T, double dt){
		
	this->T = T;
	this->N = N;
	this->rhs = rhs;
	this->dt = dt;
	
	global_t = 0.0;
	
	bodies = vector<Body>();
	toStep = vector<Body>();
		
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
		
		// then initial state
		for (int j=0; j<6; j++) {
			word = strtok(NULL, ",");
			init_state(j) = atof(word);
		}
		
		// create and store a body
		b = Body(mass, init_state);
		bodies.push_back(b);
	}
	
}
	
void NBodySolver::solve(){
		
	// while global_t < T
	// determine timesteps
	// quantize timesteps
	// determine which bodies to advance, add them to toStep
	// advanceEuler(dt, toStep)
		
}
	
void NBodySolver::advanceEuler(double dt){
	
	// create marix of values state = [[x y z vx vy vz] ... [x y z vx vy vz]]
	// take a Euler-step
	// update the bodies in toStep	
	
}

void NBodySolver::writeBodies(const char * filename){
	//iterate over bodies and write to file
	
	for (Body b: bodies){
		b.print();
	}
		
}




#include"NBodySolver.h"

	
NBodySolver::NBodySolver(int N, mat (*rhs)(vec<Body>), double T, double dt){
		
	this.T = T;
	this.N = N;
	this.rhs = rhs;
	this.dt = dt;
	
	global_t = 0.0;
	bodies = vector<Body>(1, N);
	toStep = vector<Body>();
		
}

NBodySolver::~NBodySolver(){
	
	// Destroy the system.
}
	
void NBodySolver::setInitialConditions(const char* file){
	
	Body b;
	for (int i=0; i<N; i++){
		//readline, read init pos, vel, mass
		b = Body("no_name", init_state, mass);
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
	
}




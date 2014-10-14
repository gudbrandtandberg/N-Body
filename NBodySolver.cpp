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
	
	// determine max and min timestep, and number of different
	// timesteps, based on max and min initial accelleration.
	// For example 10 different steps: [dtmax, dtmax/2, ... dtmax/2^9]
	
	while (global_t < T){

	// quantize timesteps - make lists of which bodies should have which timesteps
	// for dt in timesteps
	//   for i in numberofstepsneededtocatchup
	//	   advance(dt, toStep[dt])
		
	// also need to update times somehow...
		
		advanceRK4(dt);
		global_t += dt;
	}
	
}
	
void NBodySolver::advanceEuler(double dt){
	
	//read states from bodies
	
	for (int i=0; i<N; i++){
		states.col(i) = bodies[i].state;
	}
	
	states = states + dt*rhs(states, masses);
	
	for (int i=0; i<N; i++){
		 bodies[i].addState(states.col(i));
	}
	
}

void NBodySolver::advanceRK4(double dt){
	
	for (int i=0; i<N; i++){
		states.col(i) = bodies[i].state;
	}

	double dt2 = dt/2;
	
	mat K1 = dt*rhs(states, masses);
	mat K2 = dt*rhs(states+0.5*K1, masses);
	mat K3 = dt*rhs(states+0.5*K2, masses);
	mat K4 = dt*rhs(states+K3, masses);
	
	states = states + 1/6.0*(K1 + 2*K2 + 2*K3 + K4);
	
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




#include"NBodySolver.h"

NBodySolver::NBodySolver(int N, mat (*rhs)(mat states, vec masses), double T, double dt){
		
	this->T = T;
	this->N = N;
	this->rhs = rhs;
	this->dt = dt;
	
	global_t = 0.0;
	
	bodies = vector<Body>();
	states = zeros(6, N);
	this->masses = zeros(N);
	
	n_timesteps = 10;
	timesteps = zeros(n_timesteps);
	timesteps(0) = dt;
	for (int i=1; i<n_timesteps; i++) {
		timesteps(i) = timesteps(i-1)*2;
	}
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
		
		// create and store a body (initially with minimum timestep)
		b = Body(mass, init_state);
		b.dt = this->dt;
		bodies.push_back(b);
		
		// also store the initial states in states matrix
		states.col(i) = init_state;
	}
	
}
	
void NBodySolver::solve(){

	
	while (global_t < T){
		recomputeForces();
		advance();
		recomputeTimesteps();
		global_t += dt;

	}
	cout << timesteps << endl;
	
}

void NBodySolver::recomputeForces()
{

	for (int i=0; i<N; i++) {
		
		if ((bodies[i].nextEvalTime - global_t) < 1.E-8) {
			
			bodies[i].force = rk4(i, bodies[i].dt);;
			
		}
		else {
			//cout << "hopper over " << endl;
			continue;
		}
	
	}
}

void NBodySolver::advance(){
	
	vec nextState;
	for (int i=0; i<N; i++) {
		nextState = bodies[i].state + dt*bodies[i].force;
		bodies[i].addState(nextState);
		states.col(i) = nextState;
	}
}

void NBodySolver::recomputeTimesteps()
{
	for (int i=0; i<N; i++) {
		if ((bodies[i].nextEvalTime - global_t) < 1.E-8) {
			
			bodies[i].dt = roundBestTimestep(0.01/norm(bodies[i].force.cols(3, 5)));
			//cout << "best: " << 0.01/norm(bodies[i].force) << "rounded: " << roundBestTimestep(0.1/norm(bodies[i].force)) << endl;
			bodies[i].setNextEvalTime(global_t + bodies[i].dt);
		}
	}
}


vec NBodySolver::rk4(int i, double dt)
{
	
	mat K1 = rhs(states, masses);
	mat K2 = rhs(states + 0.5*dt*K1, masses);
	mat K3 = rhs(states + 0.5*dt*K2, masses);
	mat K4 = rhs(states + dt*K3, masses);
	
	mat forces = 1/6.0*(K1 + 2*K2 + 2*K3 + K4);
	
	return forces.col(i);
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

void NBodySolver::writeBodies(const char * filename){
	//iterate over bodies and write to file
	
	ofstream f;
	f.open(filename);
	char *buffer = new char[100];
	
	for (int body=0; body<N; body++) {
			
		f << bodies[body].state_history;
		//cout << bodies[body].state_history << endl;
		//cout << bodies[body].state_history.n_cols;
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




#define EULER 0
#define EULERCR 1
#define VERLET 2
#define RK4 3
#define AC 4   // <--- where should these be? in .h file


class ODESolver
{
	
	double T;
	int n;
	int method;
	function_pointer rhs;
	
	ODESolver((rhs(double ***), double T, int N){
		
		this.T = T;
		this.N = N;
		this.rhs = rhs;
		
	}
	
	void setInitialConditions(const char* file){
		
		//Read initial conditions from file
		
	}
	
	void solve(){
		
		// while t < T
		//	advance()
		
	}
	
	void advance(){
		
		// switch method{
		// case RK4, Euler, ...
		
	}
	
	
};
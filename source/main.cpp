/*
 * main.cpp
 *
 * Main program for Gravitational N-body simulation project by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 * Entry point for solving N-Body systems using NBodySolver.
 * First reads input parameters from the commandline, or chooses default values,
 * then initializes the solver and calls its solve method. Finally outputs the data
 * to a file.
 */

#include"NBodySolver.h"
//#include"APNBodySolver.h"

int main(int argc, char **argv)
{
	int T, N;
	int cpus = 0;
	double dt;
	double epsilon = 0.1;
	bool adaptive = 0;
	int method = 0;
	char* infile = NULL;
	
	// Read options
	if (argc > 3){
		N = atoi(argv[1]);
		T = atof(argv[2]);
		adaptive = atoi(argv[3]);
		switch (N) {
			case 3:
				infile = "./initial_conditions/solarsystem3.dat";
				break;
			case 6:
				infile = "./initial_conditions/solarsystem6.dat";
				break;
			case 11:
				infile = "./initial_conditions/solarsystem11.dat";
				break;
			case 13:
				infile = "./initial_conditions/solarsystem13.dat";
				break;
			case 18:
				infile = "./initial_conditions/solarsystem18.dat";
				break;
			case 100:
				infile = "./initial_conditions/cluster100.dat";
				break;
			case 250:
				infile = "./initial_conditions/cluster250.dat";
				break;
			case 1000:
				infile = "./initial_conditions/cluster1000.dat";
				break;
			default:
				cout << "Enter infile: ";
				cin >> infile;
				break;
		}
		
		if (adaptive) {
			if (argc > 4) {
				epsilon = atof(argv[5]);
				cpus = atoi(argv[4]);
			}
		}
		else {
			if (argc > 4){
				dt = atof(argv[4]);
				method = atoi(argv[5]);
			}
		}
		
	}
	else{
		cout << argv[0] << ": Bad usage. Should be run as either" << endl;
		cout << argv[0] << " N T adaptive" << endl;
		cout << argv[0] << " N T adaptive cpus epsilon" << endl;
		cout << argv[0] << " N T !adaptive dt method" << endl;
		exit(1);
	}
	
	// Initialize & solve
	if (!adaptive) {
		NBodySolver solver = NBodySolver(N, T, dt, method);
		solver.setInitialConditions(infile);
		solver.solve();
		solver.writeTrajectories();
		solver.writeEnergy();

	}
	else {
        /*
		APNBodySolver solver = APNBodySolver(N, T, cpus, epsilon);
		solver.setInitialConditions(infile);
		solver.solve();
		solver.writeTrajectories();
		solver.writeEnergy();*/
        int i = 0;
		
	}

}


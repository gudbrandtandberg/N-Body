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
#include"APNBodySolver.h"

int main(int argc, char **argv)
{
	int T, N;
	double dtmax, dt;
	bool adaptive = 0;
	int method = 0;
	char* infile;
	
	// Read options
	switch (argc){
	
		case 6:
			method = atoi(argv[5]);
		case 5:
			adaptive = atoi(argv[4]);
			
		case 4:
			N = atoi(argv[1]);
			T = atof(argv[2]);
			dtmax = atof(argv[3]);
			dt = atof(argv[3]);
			switch (N) {
				case 3:
					infile = "./initial_conditions/solarsystem3.dat";
					break;
				case 6:
					infile = "./initial_conditions/solarsystem6.dat";
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
				default:
					cout << "Enter infile: ";
					cin >> infile;
					break;
			}
			break;
		
		case 1:
			N = 3;
			T = 10;
			dt = 1;
			adaptive = false;
			infile = "./initial_conditions/solarsystem3.dat";
			break;
			
		default:
			cout << argv[0] << ": Bad usage. Should be run as either " << endl;
			cout << argv[0] << endl;
			cout << argv[0] << " N T dt 0 method" << endl;
			cout << argv[0] << " N T dtmax 1" << endl;
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
		APNBodySolver solver = APNBodySolver(N, T, dtmax);
		solver.setInitialConditions(infile);
		solver.solve();
		solver.writeTrajectories();
		solver.writeEnergy();
		
	}

}


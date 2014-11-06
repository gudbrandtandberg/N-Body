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

int main(int argc, char **argv)
{
	int T, N;
	double dtmax;
	bool adaptive = 0;
	char* infile;
	
	// Read options
	switch (argc){
	
		case 5:
			adaptive = atoi(argv[4]);
			
		case 4:
			N = atoi(argv[1]);
			T = atof(argv[2]);
			dtmax = atof(argv[3]);
			
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
			dtmax = 1;
			adaptive = false;
			infile = "./initial_conditions/solarsystem3.dat";
			break;
			
		default:
			cout << argv[0] << ": Bad usage. Should be run as either " << endl;
			cout << argv[0] << endl;
			cout << argv[0] << " N T dtmax" << endl;
			cout << argv[0] << " N T dtmax adaptive" << endl;
			exit(1);
	}
	
	
	// Initialize & solve
	NBodySolver solver = NBodySolver(N, T, dtmax, adaptive);
	solver.setInitialConditions(infile);
	solver.solve();
	solver.writeTrajectories();
	solver.writeEnergy();
}


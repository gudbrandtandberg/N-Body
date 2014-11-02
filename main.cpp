#include"NBodySolver.h"
#include<cmath>

using arma::vec;
using arma::mat;
using arma::norm;
using arma::zeros;
using arma::ones;

/*
 * Main program for Gravitational N-body simulation project by Gudbrand Tandberg
 * Version 1.1
 *
 * Entry point for solving N-Body systems using NBodySolver.
 * First reads input parameters from the commandline, or chooses default values,
 * then initializes the solver and calls its solve method. Finally outputs the data
 * to a file.
 */

int main(int argc, char **argv)
{
	int T, N;
	double dtmax;
	bool adaptive = 0;
	char* infile;
	char* outfile;
	
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
					infile = "./initial_conditions/solarsystem3.csv";
					outfile = "./trajectories/solarsystem3_trajectories.dat";
					break;
				case 13:
					infile = "./initial_conditions/solarsystem13.csv";
					outfile = "./trajectories/solarsystem13_trajectories.dat";
					break;
				case 18:
					infile = "./initial_conditions/solarsystem18.csv";
					outfile = "./trajectories/solarsystem18_trajectories.dat";
					break;
				case 100:
					infile = "./initial_conditions/cluster100.csv";
					outfile = "./trajectories/cluster100_trajectories.dat";
					break;
				default:
					cout << "Enter infile: ";
					cin >> infile;
					cout << "Enter outfile: ";
					cin >> outfile;
					break;
			}
			break;
		
		case 1:
			N = 3;
			T = 10;
			dtmax = 1;
			adaptive = false;
			infile = "./initial_conditions/solarsystem3.csv";
			outfile = "./trajectories/solarsystem3_trajectories.csv";
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
	solver.writeBodies(outfile);
	


}


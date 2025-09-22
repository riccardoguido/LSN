// EXERCISE 05

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"

using namespace std;

// Metropolis algorithm parameters
bool restart;						// If true restart from last position
bool state;							// Wave function state (true: Ground state, false: Excited state)
bool probability;					// Transition probability (true: Uniform, false: Gaussian)
double x0;							// Initial x position
double delta;						// Step size
int N, L;							// Number of blocks and throws in each block 
double y[3], x[3], step[3];			// Actual and proposed position
int attempted = 0, accepted = 0;	// Number of attempted and accepted moves
Random rnd;							// Random number generator

// Function to compute the statistical error
double error(double avg, double avg2, int n) {
    if (n == 0) return 0;
    return sqrt((avg2 - avg*avg)/n);
}

// Function to compute the probability density at position v
double psi2(double v[3], bool state) {
    double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	double cos_theta = v[2]/r;
	double psi;
    if (state) 
		psi = exp(-r) / sqrt(M_PI);									// Ground state (1,0,0)
	else 
		psi = 1./8. * sqrt(2./M_PI) * r * exp(-r/2.) * cos_theta;	// Excited state (2,1,0)
    return psi*psi;
}

// Function to generate a random step of the transition probability
void move(bool probability) {
    for (int i=0; i<3; i++)
        if (probability) 
			step[i] = rnd.Rannyu(-delta, delta);	// Uniform step
		else 
			step[i] = rnd.Gauss(0, delta);			// Gaussian step
}

// Function to perform one Metropolis move
void metropolis(bool state, bool probability) {
    move(probability);

    for (int i=0; i<3; i++)
        x[i] = y[i] + step[i];					 // Propose a new position

    double q = psi2(x, state) / psi2(y, state);
    double A = min(1.0, q);						 // Metropolis acceptance probability
    attempted++;								 // Count attempted move

    if (rnd.Rannyu() < A) {						 // Accept the move with probability A
        for (int j=0; j<3; j++)
            y[j] = x[j];						 // Update actual position
        accepted++;								 // Count accepted move
    }
}

int main(int argc, char *argv[]) {
	
	// Input simulation parameters reading
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) cerr << "PROBLEM: Unable to open input.dat" << endl; 

	input >> restart >> state >> probability >> x0 >> delta >> N >> L; 
    input.close();

	// Set initial position
	if (!restart) {		// Read initial position from input.dat
		y[0] = x0;
		y[1] = 0.0;
		y[2] = 0.0;
	} else {			// Read initial position from position.dat
		ifstream in_position("../INPUT/CONFIG/position.dat");
		if (!in_position.is_open()) cerr << "PROBLEM: Unable to open position.dat" << endl;

		// Read the last positions of the file
		string line, last_line;
		while (getline(in_position, line)) {
			if (line[0] != '#') last_line = line;
		}
		in_position.close();

		// Extract x, y, z from last line
		istringstream pos(last_line);
		pos >> y[0] >> y[1] >> y[2];
	}		

	// Variables for data blocking
	double avg[N] = {0.};
	double radius = 0., radius2 = 0.;
	
	// Initialize the random number generator
    rnd.SetSeed();

	// Output files
    ofstream out_radius("../OUTPUT/radius.dat");
	if(!out_radius.is_open()) cerr << "PROBLEM: Unable to open radius.dat" << endl;
	ofstream out_position("../OUTPUT/position.dat");
	if(!out_position.is_open()) cerr << "PROBLEM: Unable to open position.dat" << endl;
	ofstream out_output("../OUTPUT/output.dat");
	if(!out_output.is_open()) cerr << "PROBLEM: Unable to open output.dat" << endl;
	
	out_radius << "#      BLOCK:  	 ACTUAL_R:       R_AVE:       ERROR:" << endl;
	out_position << "#      	   X:		    Y:		  	 Z:" << endl;

    // Compute averages in each block
    for (int i=0; i<N; i++) {
		double sum = 0.;

        for (int j=0; j<L; j++) {
            metropolis(state, probability); 
			
			sum += sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);	// Distance from origin
        }
		
		// Block averages
        avg[i] = sum/L;				

		// Progressive averages and errors
		radius += avg[i];			
        radius2 += avg[i]*avg[i];

        out_radius << setw(13) << i+1	
				   << setw(13) << avg[i] 
				   << setw(13) << radius/(i+1) 
				   << setw(13) << error(radius/(i+1), radius2/(i+1), i) << endl;
				   
		out_position << setw(13) << y[0]
				     << setw(13) << y[1]
				     << setw(13) << y[2] << endl;
    }
	
	out_output << "METROPOLIS PARAMETERS:" << endl
			   << "restart = " << restart << endl
			   << "state = " << state << endl 
			   << "probability = " << probability << endl
			   << "x0 = " << x0 << endl
			   << "delta = " << delta << endl
			   << "N = " << N << endl
			   << "L = " << L << "\n" << endl
			   << "acceptance rate = " << double(accepted)/attempted << endl;
			   
    out_radius.close();
	out_position.close();
	out_output.close();
	
	// Print acceptance rate
	cout << "Acceptance rate = " << double(accepted)/attempted << endl;
	
	// Save the state of the random number generator
    rnd.SaveSeed();

    return 0;
}

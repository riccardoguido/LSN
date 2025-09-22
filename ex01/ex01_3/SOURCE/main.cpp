// EXERCISE 01_3

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"

using namespace std;

// Function to compute the statistical error
double error(double avg, double avg2, int n) {
    if (n == 0) return 0;
    return sqrt((avg2 - avg * avg) / n);
}

int main(int argc, char *argv[]) {
	
	int M = 100000;		// Number of throws
	int N = 100;		// Number of blocks
	int N_thr = M/N;	// Number of throws in each block
	
    double d = 1.0;		// Distance between the lines
    double l = 0.8;		// Length of the needle
	
	// Variables for data blocking
	double avg_pi[N] = {0.};
	double pi = 0., pi2 = 0.;
	
	// Initialize the random number generator
	Random rnd;
	rnd.SetSeed();
	
	// Output files
	ofstream out_pi("../OUTPUT/pi.dat");
	if(!out_pi.is_open()) cerr << "PROBLEM: Unable to open pi.dat" << endl;
	
	out_pi << "#      BLOCK:   ACTUAL_PI:      PI_AVE:       ERROR:" << endl;

	// Compute averages in each block
    for (int i=0; i<N; i++) {
        int N_hit = 0;
        
        for (int j=0; j<N_thr; j++) {
			// Random position of the needle center between 0 and d
            double pos = rnd.Rannyu(0., d);

			// Generate a random direction in the unitary disk
            double x, y, norm;
            do {
                x = rnd.Rannyu(-1., +1.);
                y = rnd.Rannyu(-1., +1.);
                norm = x*x + y*y;
            } while (norm > 1.);

			// Compute the cos of the angle with the x-axis
            double cos = abs(x/sqrt(norm));

			// Check if the needle crosses a line
            if (pos + 0.5*l*cos >= d or pos - 0.5*l*cos <= 0) {
                N_hit++;
            }
        }
		
		// Block averages
        avg_pi[i] = (2.*l*N_thr) / (double(N_hit)*d);

		// Progressive averages and errors
        pi += avg_pi[i];
        pi2 += avg_pi[i]*avg_pi[i];

        out_pi << setw(13) << (i+1) 
			   << setw(13) << avg_pi[i] 
			   << setw(13) << pi/(i+1)
			   << setw(13) << error(pi/(i+1), pi2/(i+1), i) << endl;
    }

    out_pi.close();
	
	// Save the state of the random number generator
	rnd.SaveSeed();
    
    return 0;
}
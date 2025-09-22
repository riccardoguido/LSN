// EXERCISE 02_2

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
    
    int M = 10000;    	// Number of random walks
    int N = 100;      	// Number of blocks
    int L = M/N;    	// Number of random walks in each block
    int steps = 100;  	// Steps per random walk
    
	// Arrays to store block averages
    double avg_discr[steps] = {0.}, avg2_discr[steps] = {0.}; 
    double avg_cont[steps] = {0.}, avg2_cont[steps] = {0.}; 
	
	// Initialize the random number generator
    Random rnd;
    rnd.SetSeed();
	
	// Output files
    ofstream out_rw_discrete("../OUTPUT/rw_discrete.dat");
    if(!out_rw_discrete.is_open()) cerr << "PROBLEM: Unable to open rw_discrete.dat" << endl;
	ofstream out_rw_continuous("../OUTPUT/rw_continuous.dat");
    if(!out_rw_continuous.is_open()) cerr << "PROBLEM: Unable to open rw_continuous.dat" << endl;
	
	out_rw_discrete << "#       STEP:  		R_AVE:       ERROR:" << endl;
	out_rw_continuous << "#       STEP:  		R_AVE:       ERROR:" << endl;

	// Compute averages in each block
    for (int i=0; i<N; i++) {
        double sum_discr[steps] = {0.};
		double sum_cont[steps] = {0.};

        for (int j=0; j<L; j++) {
            // Start from the origin
            double pos_discr[3] = {0.};	// Position for discrete random walk
			double pos_cont[3] = {0.};	// Position for continuous random walk

            // Loop over steps
            for (int k=0; k<steps; k++) {
				
                // Discrete random walk
                int dir = int(rnd.Rannyu(0.,3.));	// Direction: 0=x, 1=y, 2=z
                int sign;							// Sign: +1 or -1
				if (rnd.Rannyu() < 0.5) {
					sign = 1;
				} else {
					sign = -1;
				}
				
                pos_discr[dir] += sign;
			
                double r2_discr = pos_discr[0]*pos_discr[0] + pos_discr[1]*pos_discr[1] + pos_discr[2]*pos_discr[2]; // Distance from origin at step k
                sum_discr[k] += r2_discr;
				
				// Continuous random walk 
                double cos_theta = rnd.Rannyu(-1.,1.);				// cos_theta in [-1,1]
                double sin_theta = sqrt(1. - cos_theta*cos_theta);	// sin_theta	
                double phi = rnd.Rannyu(0.,2.*M_PI);				// phi in [0, 2*pi]

                pos_cont[0] += sin_theta * cos(phi);
                pos_cont[1] += sin_theta * sin(phi);
                pos_cont[2] += cos_theta;

                double r2_cont = pos_cont[0]*pos_cont[0] + pos_cont[1]*pos_cont[1] + pos_cont[2]*pos_cont[2]; // Distance from origin at step k
                sum_cont[k] += r2_cont;
            }
        }

        // Block averages for each step
        for (int k=0; k<steps; k++) {
            double avg_discr_step = sum_discr[k]/L;
            avg_discr[k] += avg_discr_step;
            avg2_discr[k] += avg_discr_step*avg_discr_step;
			
			double avg_cont_step = sum_cont[k]/L;
            avg_cont[k] += avg_cont_step;
            avg2_cont[k] += avg_cont_step*avg_cont_step;
        }
    }

    // Cumulative averages and errors for each step
    for (int k=0; k<steps; k++) {

        double r2_avg_discr = avg_discr[k]/N;
        double r2_avg2_discr = avg2_discr[k]/N;
        double r_discr = sqrt(r2_avg_discr);
        double err_discr = error(r2_avg_discr, r2_avg2_discr, N-1)/(2.*r_discr);
		
		double r2_avg_cont = avg_cont[k]/N;
        double r2_avg2_cont = avg2_cont[k]/N;
        double r_cont = sqrt(r2_avg_cont);
        double err_cont = error(r2_avg_cont, r2_avg2_cont, N-1)/(2.*r_cont);

        out_rw_discrete << setw(13) << k+1 
						<< setw(13) << r_discr 
						<< setw(13) << err_discr << endl;
		
        out_rw_continuous << setw(13) << k+1 
						  << setw(13) << r_cont 
						  << setw(13) << err_cont << endl;
    }

    out_rw_discrete.close();
	out_rw_continuous.close();
	
	// Save the state of the random number generator
	rnd.SaveSeed();

    return 0;
}
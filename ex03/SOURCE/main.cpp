// EXERCISE 03

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
	
    int M = 10000;     	// Number of throws
    int N = 100;       	// Number of blocks
    int L = M/N;      	// Number of throws in each block  

	// Parameters
    double S0 = 100;
    double K = 100;
	double T = 1.0;
    double r = 0.1;
    double sigma = 0.25;
	
	// Variables for data blocking
	double avg_call[N] = {0.};
	double call = 0., call2 = 0.;
	double avg_put[N] = {0.};
	double put = 0., put2 = 0.;
	
	double avg_call_d[N] = {0.};
	double call_d = 0., call2_d = 0.;
	double avg_put_d[N] = {0.};
	double put_d = 0., put2_d = 0.;

	// Initialize the random number generator
    Random rnd;
    rnd.SetSeed();

	// Output files
    ofstream out_direct("../OUTPUT/direct.dat");
	if(!out_direct.is_open()) cerr << "PROBLEM: Unable to open direct.dat" << endl;
    ofstream out_discretized("../OUTPUT/discretized.dat");
	if(!out_discretized.is_open()) cerr << "PROBLEM: Unable to open discretized.dat" << endl;
	
	out_direct << "#      BLOCK: ACTUAL_CALL:    CALL_AVE:       ERROR:  ACTUAL_PUT:     PUT_AVE:       ERROR:" << endl;
	out_discretized << "#      BLOCK: ACTUAL_CALL:    CALL_AVE:       ERROR:  ACTUAL_PUT:     PUT_AVE:       ERROR:" << endl;
	
	// Compute averages in each block
    for (int i=0; i<N; i++) {
        double sum_call = 0., sum_put = 0.;
        double sum_call_d = 0., sum_put_d = 0.;

        for (int j=0; j<L; j++) {
            // Use direct method to sample the final asset price ST
            double z = rnd.Gauss(0., 1.);
            double ST = S0*exp((r-0.5*sigma*sigma)*T + sigma*z*sqrt(T));
            double c = exp(-r*T) * max(0., ST-K);	// Direct call-option prices
            double p = exp(-r*T) * max(0., K-ST);	// Direct put-option prices

            sum_call += c;
            sum_put += p;

            // Use discretized method to sample the path of the asset price S
			double S = S0;
			int t_steps = 100;  // Time intervals
            double dt = T/t_steps;
            for (int k=0; k<t_steps; k++) {
                double z_step = rnd.Gauss(0., 1.);
                S *= exp((r-0.5*sigma*sigma)*dt + sigma*z_step*sqrt(dt));
            }
            double c_d = exp(-r*T) * max(0., S-K);	// Discretized call-option prices
            double p_d = exp(-r*T) * max(0., K-S);	// Discretized put-option prices

            sum_call_d += c_d;
            sum_put_d += p_d;
        }

        // Block averages
        avg_call[i] = sum_call/L;
        avg_put[i] = sum_put/L;
        
        avg_call_d[i] = sum_call_d/L;
        avg_put_d[i] = sum_put_d/L;
		
		// Progressive averages and errors
		call += avg_call[i];
        call2 += avg_call[i]*avg_call[i];
		put += avg_put[i];
		put2 += avg_put[i]*avg_put[i];
		
		call_d += avg_call_d[i];
        call2_d += avg_call_d[i]*avg_call_d[i];
		put_d += avg_put_d[i];
		put2_d += avg_put_d[i]*avg_put_d[i];

        out_direct << setw(13) << (i+1)
				   << setw(13) << avg_call[i]
				   << setw(13) << call/(i+1)
				   << setw(13) << error(call/(i+1), call2/(i+1), i)
				   << setw(13) << avg_put[i]
				   << setw(13) << put/(i+1)
				   << setw(13) << error(put/(i+1), put2/(i+1), i) << endl;
		
		out_discretized << setw(13) << (i+1)
					    << setw(13) << avg_call_d[i]
					    << setw(13) << call_d/(i+1)
					    << setw(13) << error(call_d/(i+1), call2_d/(i+1), i)
					    << setw(13) << avg_put_d[i]
					    << setw(13) << put_d/(i+1)
					    << setw(13) << error(put_d/(i+1), put2_d/(i+1), i) << endl;
	}
    
    out_direct.close();
    out_discretized.close();

	// Save the state of the random number generator
    rnd.SaveSeed();

    return 0;
}


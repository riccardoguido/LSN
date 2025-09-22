// EXERCISE 01_1

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
	int L = M/N;		// Number of throws in each block
	
	// Variables for data blocking
	double avg_mean[N] = {0.};
	double mean = 0., mean2 = 0.;
	
	double avg_var[N] = {0.};
	double var = 0., var2 = 0.;
	
	// Initialize the random number generator
	Random rnd;
	rnd.SetSeed();

	// Output files
	ofstream out_mean("../OUTPUT/mean.dat");
	if(!out_mean.is_open()) cerr << "PROBLEM: Unable to open mean.dat" << endl;
	ofstream out_variance("../OUTPUT/variance.dat");
	if(!out_variance.is_open()) cerr << "PROBLEM: Unable to open variance.dat" << endl;
	
	out_mean << "#      BLOCK: ACTUAL_MEAN:    MEAN_AVE:       ERROR:" << endl;
	out_variance << "#      BLOCK:  ACTUAL_VAR:     VAR_AVE:       ERROR:" << endl;
	
	// Compute averages in each block
    for (int i=0; i<N; i++) {
        double sum_mean = 0.;
        double sum_var = 0.;
		
        for (int j=0; j<L; j++) {
            double r = rnd.Rannyu();
            sum_mean += r;
            sum_var += (r - 0.5)*(r - 0.5);
        }
		
		// Block averages
        avg_mean[i] = sum_mean / L;			
        avg_var[i] = sum_var / L;
		
        // Progressive averages and errors
		mean += avg_mean[i];
        mean2 += avg_mean[i]*avg_mean[i];
		
        var += avg_var[i];
        var2 += avg_var[i]*avg_var[i];			

        out_mean << setw(13) << (i+1) 
				 << setw(13) << avg_mean[i] 
				 << setw(13) << mean/(i+1)
				 << setw(13) << error(mean/(i+1), mean2/(i+1), i) << endl;
				 
        out_variance << setw(13) << (i+1) 
					 << setw(13) << avg_var[i] 
					 << setw(13) << var/(i+1) 
					 << setw(13) << error(var/(i+1), var2/(i+1), i) << endl;
    }

    out_mean.close();
	out_variance.close();

// ==================== CHI^2 TEST ====================
	
    int m = 100;    		// Number of sub-intervals
    int n = 10000; 			// Number of throws per experiment
    int N_tests = 100;  	// Number of chi^2 repetitions
	
	// Output file
	ofstream out_chi("../OUTPUT/chi.dat");
	if(!out_chi.is_open()) cerr << "PROBLEM: Unable to open chi.dat" << endl;
	
	out_chi << "#       TEST:  	    CHI^2:" << endl;

	// Loop over N_tests experiments
    for (int i=0; i<N_tests; i++) {
        double chi = 0.;
        int n_i[m] = {}; 

		// Generate uniform random numbers
        for (int j=0; j<n; j++) {
            double r = rnd.Rannyu();
			int bin = int(r*m);
            n_i[bin]++;
        }

        // Compute chi^2
		double np = n/m;
        for (int k=0; k<m; k++) {
            chi += (n_i[k] - np)*(n_i[k] - np) / np;
        }

        out_chi << setw(13) << (i+1) 
				<< setw(13) << chi << endl;
    }

    out_chi.close();
	
	// Save the state of the random number generator
	rnd.SaveSeed();
	
    return 0;
}
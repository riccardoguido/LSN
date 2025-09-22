// EXERCISE 02_1

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
	
    int M = 1000000;	// Number of throws
	int N = 100;		// Number of blocks
	int L = M/N;		// Number of throws in each block
	
	// Variables for data blocking
	double avg[N] = {0.};
	double integral = 0., integral2 = 0.;
	
	double avg_is[N] = {0.};
	double integral_is = 0., integral2_is = 0.;
	
	// Initialize the random number generator
	Random rnd;
	rnd.SetSeed();
	
	// Output files
	ofstream out_integral("../OUTPUT/integral.dat");
	if(!out_integral.is_open()) cerr << "PROBLEM: Unable to open integral.dat" << endl;
	ofstream out_integral_is("../OUTPUT/integral_is.dat");
	if(!out_integral_is.is_open()) cerr << "PROBLEM: Unable to open integral_is.dat" << endl;
	
	out_integral << "#      BLOCK:    ACTUAL_I:       I_AVE:       ERROR:" << endl;
	out_integral_is << "#      BLOCK:    ACTUAL_I:       I_AVE:       ERROR:" << endl;

    // Compute averages in each block
    for (int i=0; i<N; i++) {
        double sum = 0.;
		double sum_is = 0.;
		
        for (int j=0; j<L; j++) {
			// Use uniforme distribution of x to estimate the integral
            double x = rnd.Rannyu();
			double f = 0.5*M_PI*cos(0.5*M_PI*x);
            sum += f;
			// Use importance sampling of x with a linear distribution p(x) = 2(1 - x) to estimate the integral
			double x_is = rnd.Linear();
			double f_is = 0.5*M_PI*cos(0.5*M_PI*x_is); 
			double p = 2.*(1.-x_is);	// PDF evaluated at x_is
            sum_is += f_is/p;
        }
		
		// Block averages
        avg[i] = sum/L;						
        avg_is[i] = sum_is/L;		
		
		// Progressive averages and errors
		integral += avg[i];
        integral2 += avg[i]*avg[i];
		
		integral_is += avg_is[i];
        integral2_is += avg_is[i]*avg_is[i];	
		
		out_integral << setw(13) << (i+1)
					 << setw(13) << avg[i]
					 << setw(13) << integral/(i+1) 
					 << setw(13) << error(integral/(i+1), integral2/(i+1), i) << endl;
					 
		out_integral_is << setw(13) << (i+1)
					    << setw(13) << avg_is[i]
						<< setw(13) << integral_is/(i+1) 
						<< setw(13) << error(integral_is/(i+1), integral2_is/(i+1), i) << endl;
    }
	
	out_integral.close();
	out_integral_is.close();
	
	// Save the state of the random number generator
	rnd.SaveSeed();

    return 0;
}
// EXERCISE 01_2

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
	
	int M = 10000; 				// Number of throws
	int N[4] = {1, 2, 10, 100};	// Number of experiment repetitions
	
	// Initialize the random number generator
	Random rnd;
	rnd.SetSeed();
	
	// Output files
	ofstream out_uniform("../OUTPUT/uniform.dat");
	if(!out_uniform.is_open()) cerr << "PROBLEM: Unable to open uniform.dat" << endl;
	ofstream out_exponential("../OUTPUT/exponential.dat");
	if(!out_exponential.is_open()) cerr << "PROBLEM: Unable to open exponential.dat" << endl;
	ofstream out_lorentzian("../OUTPUT/lorentzian.dat");
	if(!out_lorentzian.is_open()) cerr << "PROBLEM: Unable to open lorentzian.dat" << endl;
	
	out_uniform << "#        S_1:  	      S_2:        S_10:       S_100:" << endl;
	out_exponential << "#        S_1:  	      S_2:        S_10:       S_100:" << endl;
	out_lorentzian << "#        S_1:  	      S_2:        S_10:       S_100:" << endl;
	
	// Loop over M experiments
	for(int i=0; i<M; i++) {
		double unif[4] = {0.};		// Uniform distribution
		double exp[4] = {0.};		// Exponential distribution
		double lorentz[4] = {0.};	// Lorentzian distribution
		
		// Compute sums of N samples
		for (int j=0; j<4; j++) {
			
			for(int k=0; k<N[j]; k++) {
				unif[j] += rnd.Rannyu();
				exp[j] += rnd.Exp(1);
				lorentz[j] += rnd.Lorentz(0,1);
			}
		}
		
		// Averges of N samples
		out_uniform << setw(13) << unif[0]/N[0] 
					<< setw(13) << unif[1]/N[1] 
					<< setw(13) << unif[2]/N[2] 
					<< setw(13) << unif[3]/N[3] << endl;
				 
		out_exponential << setw(13) << exp[0]/N[0] 
						<< setw(13) << exp[1]/N[1] 
						<< setw(13) << exp[2]/N[2] 
						<< setw(13) << exp[3]/N[3] << endl;
		
		out_lorentzian << setw(13) << lorentz[0]/N[0] 
					   << setw(13) << lorentz[1]/N[1] 
					   << setw(13) << lorentz[2]/N[2] 
					   << setw(13) << lorentz[3]/N[3] << endl;
	}
	
	out_uniform.close();
	out_exponential.close();
	out_lorentzian.close();
	
	// Save the state of the random number generator
	rnd.SaveSeed();
	
	return 0;
}
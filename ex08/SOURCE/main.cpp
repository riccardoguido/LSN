// EXERCISE 08

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
double mu, sigma;					// Wave function parameters 
double x0;							// Initial x position
double delta;						// Step size
int N, L;							// Number of blocks and throws in each block 
double y, x;						// Actual and proposed position
int attempted = 0, accepted = 0;	// Number of attempted and accepted moves

// Simulated Annealing parameters
bool annealing;						// If true do Simulated Annealing
double beta, dbeta;					// Inverse temperature and its increment
int steps_beta;						// Number of optimization steps for each beta value
double delta_mu, delta_sigma;		// Maximum variation of mu and sigma during a step
double err_mu, err_sigma;			// Target error for mu and sigma

Random rnd;							// Random number generator	

// Function to compute the statistical error
double error(double avg, double avg2, int n) {
    if (n == 0) return 0;
    return sqrt((avg2 - avg*avg)/n);
}

// Function to compute the probability density at position v
double psi2(double v) {
    double exp1 = exp(-pow(v-mu,2) / (2*sigma*sigma));
    double exp2 = exp(-pow(v+mu,2) / (2*sigma*sigma));
	double psi = exp1 + exp2;
    return psi*psi;
}

// Function to compute the local energy at position v
double energy_loc(double v) {
    double exp1 = exp(-pow(v-mu,2) / (2*sigma*sigma));
    double exp2 = exp(-pow(v+mu,2) / (2*sigma*sigma));
	double psi = exp1 + exp2;
    double K = -0.5 * (exp1 * (pow(v-mu,2) / pow(sigma,4) - 1.0 / (sigma*sigma))
               + exp2 * (pow(v+mu,2) / pow(sigma,4) - 1.0 / (sigma*sigma)));
    double V = pow(v,4) - 2.5*pow(v,2);
    return K/psi + V;
}

// Function to perform one Metropolis move
void metropolis() {
	double step = rnd.Rannyu(-delta, delta); // Uniform step
    x = y + step;							 // Propose a new position
	
    double q = psi2(x) / psi2(y);			
    double A = min(1.0, q);					 // Metropolis acceptance probability
    attempted++;							 // Count attempted move

    if (rnd.Rannyu() < A) {					 // Accept the move with probability A
        y = x;								 // Update actual position
        accepted++;							 // Count accepted move
    }
}

// Function to compute the expectation value for Hamiltonian
double hamiltonian(double& err, bool print_blocks = false) {
	if (!restart) {		// Read initial position from input.dat
		y = x0;
	} else {			// Read initial position from position.dat
		ifstream in_position("../INPUT/CONFIG/position.dat");
		if (!in_position.is_open()) cerr << "PROBLEM: Unable to open position.dat" << endl;
		string line, last_line;
		while (getline(in_position, line)) {
			if (line[0] != '#') last_line = line;
		}
		in_position.close();
		istringstream pos(last_line);
		pos >> y;
	}

	double H_block = 0, H = 0., H2 = 0.;	 // Variables for data blocking
	
	ofstream out_energy, out_position;		 // Output files
	if (print_blocks) {						 
		out_energy.open("../OUTPUT/energy.dat");
		if (!out_energy.is_open()) cerr << "PROBLEM: Unable to open energy.dat" << endl;
		out_position.open("../OUTPUT/position.dat");
		if(!out_position.is_open()) cerr << "PROBLEM: Unable to open position.dat" << endl;
		out_energy << "#      BLOCK:  	 ACTUAL_H:       H_AVE:       ERROR:" << endl;
		out_position << "#      	   X:" << endl;
	}

    for (int i=0; i<N; i++) {				 // Compute averages in each block
		double sum = 0.;
        for (int j=0; j<L; j++) {
            metropolis(); 	
			sum += energy_loc(y);
        }
		H_block = sum/L;				
		H += H_block;			
        H2 += H_block*H_block;
		
		if (print_blocks) {
			out_energy << setw(13) << i+1
					   << setw(13) << H_block 
					   << setw(13) << H/(i+1) 
					   << setw(13) << error(H/(i+1), H2/(i+1), i) << endl;	   
			out_position << setw(13) << y << endl;	 
		}
    }
	if (print_blocks) {
		out_energy.close();
		out_position.close();
	}
	
    H /= N;									 // Cumulative average
    H2 /= N;
    err = error(H, H2, N-1);				 // Cumulative error
    return H;
}

int main(int argc, char *argv[]) {
	
	// Input simulation parameters reading
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) cerr << "PROBLEM: Unable to open input.dat" << endl; 

	input >> restart >> mu >> sigma >> x0 >> delta >> N >> L 
		  >> annealing >> beta >> dbeta >> steps_beta >> delta_mu >> delta_sigma >> err_mu >> err_sigma;
    input.close();
	
	// Initialize the random number generator
    rnd.SetSeed();
	
	// Output file
	ofstream out_output("../OUTPUT/output.dat");
	if(!out_output.is_open()) cerr << "PROBLEM: Unable to open output.dat" << endl;
	
	double err;
	
	if (!annealing) {
		
		hamiltonian(err, true);
		
		out_output << "METROPOLIS PARAMETERS:" << endl 
				   << "restart = " << restart << endl 
				   << "mu = " << mu << endl 
				   << "sigma = " << sigma << endl
				   << "x0 = " << x0 << endl
				   << "delta = " << delta << endl
				   << "N = " << N << endl
				   << "L = " << L << "\n" << endl
				   << "acceptance rate = " << double(accepted)/attempted << endl;
			   
		out_output.close();
		
		// Print acceptance rate
		cout << "Acceptance rate = " << double(accepted)/attempted << endl;
	}

// =============== SIMULATED ANNEALING ALGORITHM: variational optimization of energy H ===============
	if (annealing) {
		
		// Output file
        ofstream out_annealing("../OUTPUT/annealing.dat");		
        if (!out_annealing.is_open()) cerr << "PROBLEM: Unable to open annealing.dat" << endl;
		
		out_annealing << "#       STEP:  	 	 BETA:	  		MU:       SIGMA:       H_AVE:       ERROR:" << endl;
		out_output << "METROPOLIS PARAMETERS:" << endl 
				   << "restart = " << restart << endl 
				   << "mu = " << mu << endl 
				   << "sigma = " << sigma << endl
				   << "x0 = " << x0 << endl
				   << "delta = " << delta << endl
				   << "N = " << N << endl
				   << "L = " << L << "\n" << endl
				   << "SIMULATED ANNEALING PARAMETERS:" << endl
				   << "annealing = " << annealing << endl 
				   << "beta = " << beta << endl   
				   << "dbeta = " << dbeta << endl      
				   << "steps_beta = " << steps_beta << endl   
				   << "delta_mu = " << delta_mu << endl	  
				   << "delta_sigma = " << delta_sigma << endl
				   << "err_mu = " << err_mu << endl		  
				   << "err_sigma = " << err_sigma << "\n" << endl;

		// Initial energy estimate
		double H_old = hamiltonian(err);
		
		// Progressively increase beta to cool the system
        while (2*delta_mu/beta > err_mu or 2*delta_sigma/beta > err_sigma) {
			
            for (int i=0; i<steps_beta; i++) {
				// Propose small changes of mu and sigma
                double dmu = rnd.Rannyu(-delta_mu, delta_mu) / beta;
                double dsigma = rnd.Rannyu(-delta_sigma, delta_sigma) / beta;	
                double mu_new = mu + dmu;
                double sigma_new = sigma + dsigma;
				
				// Backup current parameters in case of rejection
                double mu_backup = mu;
                double sigma_backup = sigma;
				
				// Compute energy in the new state
                mu = mu_new;
                sigma = sigma_new;
                double H_new = hamiltonian(err);
				
				// Energy difference
                double dH = H_old - H_new;			
                double A = min(1.0, exp(beta*dH));
                if (rnd.Rannyu() < A) {		// Accept new parameters
                    H_old = H_new;			
                } else {					// Reject move and restore previous state
                    mu = mu_backup;
                    sigma = sigma_backup;	
                }

                out_annealing << setw(13) << i+1
							  << setw(13) << beta	
							  << setw(13) << mu 
							  << setw(13) << sigma 
							  << setw(13) << H_old 
							  << setw(13) << err << endl;
            }
			
			// Cool down the system
            beta += dbeta;
        }
		
		out_output << "OPTIMIZATION RESULTS:" << endl 
				   << "H = " << H_old << endl
				   << "mu = " << mu << endl 
				   << "sigma = " << sigma << endl;
				   
		out_annealing.close();
		out_output.close();
		
        cout << "H = " << H_old << endl;
        cout << "mu = " << mu << "  sigma = " << sigma << endl;
	}
			   
	// Save the state of the random number generator
    rnd.SaveSeed();

    return 0;
}
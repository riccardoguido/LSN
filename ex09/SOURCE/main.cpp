// EXERCISE 09

#include <fstream>
#include <iomanip>
#include "genetic.h"

using namespace std;

int main(int argc, char *argv[]) {
	
	// Problem parameters
	int Ncities;   	 // Number of cities
    int Npop;   	 // Population size
    int Ngen;  		 // Number of generations

    // Genetic algorithm parameters
	bool type;		 // Cities generation type (true: Circle, false: Square)
    double psel;	 // Rank selection exponent
    double pcross;   // Crossover probability
    double pm_swap;  // Swap mutation probability
    double pm_inv;   // Inversion mutation probability
    double pm_shift; // Shift mutation probability
    double pm_mperm; // Permutation mutation probability
	
    // Input simulation parameters reading
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) cerr << "PROBLEM: Unable to open input.dat" << endl; 

	input >> Ncities >> Npop >> Ngen
		  >> type >> psel >> pcross >> pm_swap >> pm_inv >> pm_shift >> pm_mperm; 
    input.close();

    // Initialize the random number generator
    Random rnd;
    rnd.SetSeed();
	
	// Vector of city coordinates
	vector<City> cities;
	
	// Generate city positions
	if (type)
		generate_circular_cities(cities, Ncities, rnd);  // Generate Ncities cities on a circle
	if (!type)
		generate_square_cities(cities, Ncities, rnd);	 // Generate Ncities cities in a square
		
	// Initialize the population of chromosomes
	Population pop = init_population(Npop, Ncities, cities, rnd);
	// Sort population by fitness
	sort_population(pop);

	// Output files
	ofstream out_length("../OUTPUT/best_length.dat");
	if(!out_length.is_open()) cerr << "PROBLEM: Unable to open best_length.dat" << endl;
	ofstream out_half_length("../OUTPUT/avg_best_half_length.dat");
	if(!out_half_length.is_open()) cerr << "PROBLEM: Unable to open avg_best_half_length.dat" << endl;
	ofstream out_path("../OUTPUT/best_path.dat");
	if(!out_path.is_open()) cerr << "PROBLEM: Unable to open best_path.dat" << endl;
	ofstream out_output("../OUTPUT/output.dat");
	if(!out_output.is_open()) cerr << "PROBLEM: Unable to open output.dat" << endl;
	
	out_length << "#      	 GEN:      LENGTH:" << endl;
	out_half_length << "#      	 GEN:      LENGTH:" << endl;
	out_path << "#      	CITY:			X:     		 Y:" << endl;
	out_output << "PROBLEM PARAMETERS:" << endl
			   << "Ncities = " << Ncities << endl 
			   << "Npop = " << Npop << endl
			   << "Ngen = " << Ngen << "\n" << endl
			   << "GENETIC PARAMETERS:" << endl
			   << "type = " << type << endl
			   << "psel = " << psel << endl
			   << "pcross = " << pcross << endl
			   << "pm_swap = " << pm_swap << endl
			   << "pm_inv = " << pm_inv << endl
			   << "pm_shift = " << pm_shift << endl
			   << "pm_mperm = " << pm_mperm << endl;
	
	// Loop over Ngen generations
	for (int i=0; i<Ngen; i++) {
		evolve_one_generation(pop, cities, rnd, psel, pcross, pm_swap, pm_inv, pm_shift, pm_mperm);

		// Best path length
		out_length << setw(13) << (i+1)
				   << setw(13) << best_length(pop) << endl;

		// Average path length of best half population
		out_half_length << setw(13) << (i+1)
				        << setw(13) << avg_best_half_length(pop) << endl;

		// Coordinates of best path in cartesian coordinates
		const Chromosome& best = pop.individuals[0];
		
		for (int j=0; j<Ncities; j++) {
			int idx = best.path[j] - 1;  // Path is 1-based
			out_path << setw(13) << (j+1) 
					 << setw(13) << cities[idx].x 
				     << setw(13) << cities[idx].y << endl;
		}
	}
	
	out_length.close();
	out_half_length.close();
	out_path.close();
	out_output.close();
	
	// Save the state of the random number generator
    rnd.SaveSeed();

    return 0;
}

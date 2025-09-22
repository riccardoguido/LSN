// EXERCISE 10

#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <mpi.h>
#include "genetic.h"

using namespace std;

int main(int argc, char* argv[]) {
	
	// Problem parameters
	int Ncities;		// Number of cities
    int Npop;   	 	// Population size for each continent
    int Ngen;  		 	// Number of generations
    
	// Message Passing Interface (MPI) parameters
	int size;  			// Total number of processes (continents)
    int rank;  			// Process identifier
	bool migration;  	// If true do migration between continents
	int Nmigr;   	    // Number of generations between every migration
	
    // Genetic algorithm parameters
    double psel;	 	// Rank selection exponent
    double pcross;   	// Crossover probability
    double pm_swap;  	// Swap mutation probability
    double pm_inv;   	// Inversion mutation probability
    double pm_shift; 	// Shift mutation probability
    double pm_mperm; 	// Permutation mutation probability
	
    // Input simulation parameters reading
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) cerr << "PROBLEM: Unable to open input.dat" << endl; 

	input >> Npop >> Ngen 
		  >> migration >> Nmigr
		  >> psel >> pcross >> pm_swap >> pm_inv >> pm_shift >> pm_mperm;
    input.close();
	
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// Initialize the random number generator
    Random rnd;
    rnd.SetSeed();                  
    for (int i=0; i<1000*rank; i++) rnd.Rannyu(); // Each process gets a different random sequence
 
	// Vector of city coordinates
    vector<City> cities;

	// Rank 0 reads city coordinates from input file
    if (rank == 0) {
        load_cities(cities, "../INPUT/cap_prov_ita.dat");
		Ncities = (int)cities.size();
    } else {
        Ncities = 0;
    }
	
	// Rank 0 broadcasts Ncities to all other ranks
    MPI_Bcast(&Ncities, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        cities.resize(Ncities);
    }
	
    // Rank 0 broadcasts city coordinates to all other ranks
    vector<double> xs(Ncities);
    vector<double> ys(Ncities);

    if (rank == 0) {
        for (int i=0; i<Ncities; i++) {
            xs[i] = cities[i].x;
            ys[i] = cities[i].y;
        }
    }
    MPI_Bcast(xs.data(), Ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ys.data(), Ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i=0; i<Ncities; i++) {
        cities[i].x = xs[i];
        cities[i].y = ys[i];
    }

    // Initialize the population of chromosomes
	Population pop = init_population(Npop, Ncities, cities, rnd);
	// Sort population by fitness
	sort_population(pop);

    // Vectors used during migration
    vector<int> send_path(Ncities); 	// Store the best path to send  
    vector<int> recv_path(Ncities); 	// Store the best path received 
	
	// Output files
	ostringstream dirname;
	dirname << "../OUTPUT/rank_" << rank;
	string dirpath = dirname.str();
	filesystem::create_directories(dirpath);
	
	ofstream out_length(dirpath + "/best_length.dat");
	if(!out_length.is_open()) cerr << "PROBLEM: Unable to open best_length.dat" << endl;
	ofstream out_half_length(dirpath + "/avg_best_half_length.dat");
	if(!out_half_length.is_open()) cerr << "PROBLEM: Unable to open avg_best_half_length.dat" << endl;
	ofstream out_path(dirpath + "/best_path.dat");
	if(!out_path.is_open()) cerr << "PROBLEM: Unable to open best_path.dat" << endl;
	
	out_length << "#      	 GEN:      LENGTH:" << endl;
	out_half_length << "#      	 GEN:      LENGTH:" << endl;
	out_path << "#      	CITY:			X:     		 Y:" << endl;

	// Decide whether to perform a migration
    for (int i=0; i<Ngen; i++) {
        bool do_migration = false;
        if (migration) {
            if (((i+1)%Nmigr == 0) and (size > 1)) {
                do_migration = true;
            }
        }
		
		// Indipendent evolution
        if (!do_migration) evolve_one_generation(pop, cities, rnd, psel, pcross, pm_swap, pm_inv, pm_shift, pm_mperm);
		
		// Random migration between continents
        if (do_migration) {
       
			// Rank 0 selects two distinct processes for migration: giver and receiver
            int giver = 0, receiver = 0;
			
            if (rank == 0) {
                giver = (int)rnd.Rannyu(0, size);
					do {
					receiver = (int)rnd.Rannyu(0, size);
				} while (receiver == giver);
			}     
            MPI_Bcast(&giver, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&receiver, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // giver sends its best path (firt index)
            if (rank == giver) {
                send_path = pop.individuals[0].path; 
                MPI_Send(send_path.data(), Ncities, MPI_INT, receiver, 0, MPI_COMM_WORLD);
            }

            // receiver replaces its worst path (last index) with the received one
            if (rank == receiver) {
                MPI_Recv(recv_path.data(), Ncities, MPI_INT, giver, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pop.individuals[Npop - 1].path = recv_path; 
                pop.individuals[Npop - 1].check();                       // Validate the new path
				evaluate_fitness(cities, pop.individuals[Npop - 1]);     // Compute length and fitness of the new path
				sort_population(pop);                                 	 // Reorder population
            }
        }
		
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

    // Find the best global length among all ranks
    struct {double val; int rank;} in_minloc, out_minloc;
	in_minloc.val  = best_length(pop);
	in_minloc.rank = rank;
	MPI_Reduce(&in_minloc, &out_minloc, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

	// Rank 0 writes the parallelization results
	if (rank == 0) {
		ofstream out_output("../OUTPUT/output.dat");
		if(!out_output.is_open()) cerr << "PROBLEM: Unable to open output.dat" << endl;
		
		out_output << "PROBLEM PARAMETERS:" << endl
				   << "Npop = " << Npop << endl
				   << "Ngen = " << Ngen << "\n" << endl
				   << "MPI PARAMETERS:" << endl
				   << "migration = " << migration << endl
				   << "Nmigr = " << Nmigr << "\n" << endl
				   << "GENETIC PARAMETERS:" << endl
				   << "psel = " << psel << endl
				   << "pcross = " << pcross << endl
				   << "pm_swap = " << pm_swap << endl
				   << "pm_inv = " << pm_inv << endl
				   << "pm_shift = " << pm_shift << endl
				   << "pm_mperm = " << pm_mperm << "\n" << endl
				   << "PARALLELIZATION RESULTS:" << endl
				   << "total number of processes = " << size << endl
				   << "best global length = " << out_minloc.val << " found by rank " << out_minloc.rank << endl;
				   
		out_output.close();
		
		cout << "Total number of processes = " << size << endl;
		cout << "Best global length = " << out_minloc.val << " found by rank " << out_minloc.rank << endl;
	}
	
	out_length.close();
	out_half_length.close();
	out_path.close();
	
	// Save the state of the random number generator
    rnd.SaveSeed();
	
	// Finalize the MPI environment
    MPI_Finalize();
	
    return 0;
}
#ifndef __Genetic__
#define __Genetic__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include "random.h"

using namespace std;

// Structure to represent a city in cartesian coordinates
struct City {
    double x, y;
};

// Structure to represent a chromosome
struct Chromosome {
    vector<int> path;  		 // Permutation of cities with city 1 fixed at the head
    double fitness = 0.0;    // Selection score: 1/L(2)
    double length  = 0.0;    // Total path length: L(1)
    void check() const;      // Validate the chromosome: path[0]==1 and path is a true permutation
};

// Structure to represent a population of chromosomes
struct Population {
    vector<Chromosome> individuals;
};

// ==================== CITY GENERATION ====================

// Method to generate N cities uniformly on the unit circle
void generate_circular_cities(vector<City>& cities, int N, Random& rnd);
// Method to generate N cities uniformly in the square
void generate_square_cities(vector<City>& cities, int N, Random& rnd);
// Method to load city coordinates from file
void load_cities(vector<City>& cities, const string& filename);

// ==================== CHROMOSOME & FITNESS ====================

// Method to create a random chromosome of size N with city 1 fixed at the head
Chromosome make_random_chromosome(int N, Random& rnd);
// Method to compute length and fitness of a chromosome
void evaluate_fitness(vector<City>& cities, Chromosome& c);

// ==================== POPULATION ====================

// Method to initialize a population of size M with random chromosomes
Population init_population(int M, int N, vector<City>& cities, Random& rnd);
// Method to compare two chromosomes by fitness
bool compare_chromosomes(const Chromosome& a, const Chromosome& b);
// Method to sort a population from best to worst by fitness
void sort_population(Population& pop);

// ==================== SELECTION ====================

// Method to select one index by rank selection
int select_rank(Population& pop, Random& rnd, double p);
// Method to select two distinct parent indices by rank selection
void select_two_parents(Population& pop, Random& rnd, double p, int& i1, int& i2);

// ==================== MUTATIONS ====================

// Method to swap two positions in the tail
void mutate_swap (Chromosome& c, Random& rnd);
// Method to invert a contiguous segment in the tail
void mutate_inversion (Chromosome& c, Random& rnd);
// Method to shift a contiguous block in the tail to the right by n positions
void mutate_shift (Chromosome& c, Random& rnd);
// Method to swap two disjoint blocks of equal length m in the tail
void mutate_mpermut (Chromosome& c, Random& rnd);

// ==================== CROSSOVER ====================

// Method to perform order-based crossover on a segment
void crossover(Chromosome& child, Chromosome& parent2, int pos, int len);

// ==================== EVOLUTION & METRICS ====================

// Method to evolve the population by one generation
void evolve_one_generation(Population& pop, vector<City>& cities, Random& rnd,
                           double psel, double pcross, double pm_swap, double pm_inv, double pm_shift, double pm_mperm);

// Method to get the best path length in the population
double best_length(Population& pop);
// Method to get the average path length over the best half of the population
double avg_best_half_length(Population& pop);

#endif // __Genetic__

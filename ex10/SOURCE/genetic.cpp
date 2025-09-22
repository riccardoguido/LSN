#include "genetic.h"

using namespace std;

void Chromosome::check() const{
	// Check path is not empty and city 1 is fixed at the head
    if (path.empty() or path[0]!=1) {
        cerr << "PROBLEM: First city must be 1" << endl; 
		exit(1);
    }
	// Check path is a true permutation
    vector<int> tmp = path;
    sort(tmp.begin(), tmp.end());
    for (int i = 0; i < static_cast<int>(tmp.size()); i++) {
        if (tmp[i] != i+1) {
            cerr << "PROBLEM: Invalid permutation" << endl; 
			exit(2);
        }
    }
}

// ==================== CITY GENERATION ====================

void generate_circular_cities(vector<City>& cities, int N, Random& rnd) {
    cities.resize(N);
    for (int i=0; i<N; i++) {
        double phi = rnd.Rannyu(0.0, 2.0*M_PI);
        cities[i] = {cos(phi), sin(phi)};
    }
}

void generate_square_cities(vector<City>& cities, int N, Random& rnd) {
    cities.resize(N);
    for (int i=0; i<N; i++) {
        double x = rnd.Rannyu(-1.0, 1.0);
        double y = rnd.Rannyu(-1.0, 1.0);
        cities[i] = {x, y};
    }
}

void load_cities(vector<City>& cities, const string& filename) {
    ifstream input(filename);
    if (!input.is_open()) cerr << "PROBLEM: Unable to open " << filename << endl;
    cities.clear();
    City c;
    while (input >> c.x >> c.y) cities.push_back(c);
    input.close();
}

// ==================== CHROMOSOME & FITNESS ====================

Chromosome make_random_chromosome(int N, Random& rnd) {
    Chromosome c;
	// Initialize path as [1,2,...,N]
    c.path.resize(N);
    for (int i=0; i<N; i++) 
		c.path[i] = i+1;
	// Swap only indices 1..N-1 (keep city 1 fixed)
    for (int i=N-1; i>=2; i--) {
        int j = (int)rnd.Rannyu(1, i+1); 	// j in [1,i], never swaps index 0
        swap(c.path[i], c.path[j]);
    }
	// Check validity
    c.check();
    return c;
}

void evaluate_fitness(vector<City>& cities, Chromosome& c) {
    int N = c.path.size();
    double sum_len = 0.0;
    double sum_len_sq = 0.0;

    for (int i=0; i<N; i++) {
		int a = c.path[i] - 1;           // Convert 1-based index to 0-based
        int b = c.path[(i+1) % N] - 1;   // Connect last to first
		
		// Compute squared distance
        double dx = cities[b].x - cities[a].x;
        double dy = cities[b].y - cities[a].y;
        double d2 = dx*dx + dy*dy;
		// Add to totals
        sum_len += sqrt(d2);			 // L(1)
        sum_len_sq += d2;			     // L(2)
    }
	// Save total length and fitness
    c.length = sum_len;					 
    c.fitness = 1.0 / sum_len_sq;		
}

// ==================== POPULATION ====================

Population init_population(int M, int N, vector<City>& cities, Random& rnd) {
    Population pop;
    pop.individuals.resize(M);
    for (int i=0; i<M; i++) {
        pop.individuals[i] = make_random_chromosome(N, rnd);
        evaluate_fitness(cities, pop.individuals[i]);
    }
    return pop;
}

bool compare_chromosomes(const Chromosome& a, const Chromosome& b) {
    return a.fitness > b.fitness; 	
}

void sort_population(Population& pop) {
	// Sort in descending order (best individuals first)
    sort(pop.individuals.begin(), pop.individuals.end(), compare_chromosomes);
}

// ==================== SELECTION ====================

int select_rank(Population& pop, Random& rnd, double p) {
    int M = pop.individuals.size();
    double r = rnd.Rannyu();
    int idx = (int)((M - 1) * pow(r, p));
    return idx;
}

void select_two_parents(Population& pop, Random& rnd, double p, int& i1, int& i2) {
    i1 = select_rank(pop, rnd, p);
    do {
		i2 = select_rank(pop, rnd, p); 
		} while (i2 == i1);
}

// ==================== MUTATIONS ====================

void mutate_swap(Chromosome& c, Random& rnd) {
    int N = c.path.size();
    int i = 0, j = 0;
    i = 1 + (int)rnd.Rannyu(0, N-1);      // i in [1, N-1]
    do {
        j = 1 + (int)rnd.Rannyu(0, N-1);  // j in [1, N-1], distinct from i
    } while (j == i);
	
	// Swap the two cities
    swap(c.path[i], c.path[j]);
}

void mutate_inversion(Chromosome& c, Random& rnd) {
    int N = c.path.size();
    int pos = 1 + (int)rnd.Rannyu(0, N-2);    // Position in [1, N-2]
    int max_m = N - pos;                      // Maximum possible segment length
    int m = 2 + (int)rnd.Rannyu(0, max_m-1);  // Block length
    int L = pos + m - 1;					  // End index of the segment

    for (int a=pos, b=L; a<b; a++, b--)
		// Swap elements between pos and L
        swap(c.path[a], c.path[b]);
}

void mutate_shift(Chromosome& c, Random& rnd) {
    int N = c.path.size();  
    int pos = 1 + (int)rnd.Rannyu(0, N-1);	  // Position in [1, N-1]
    int m = 1 + (int)rnd.Rannyu(0, N-pos);	  // Block length
    int max_n = (N - 1) - (pos + m - 1);	  // Maximum possible shift distance
    int n;
	
    if (max_n > 0)
        n = 1 + (int)rnd.Rannyu(0, max_n);
    else 
        n = 0; 
    if (n == 0) return;
	
    // Copy block
    vector<int> block(c.path.begin() + pos, c.path.begin() + pos + m);
	// Remove block
    c.path.erase(c.path.begin() + pos, c.path.begin() + pos + m);
	// Reinsert block shifted by n
    c.path.insert(c.path.begin() + pos + n, block.begin(), block.end());
}

void mutate_mpermut(Chromosome& c, Random& rnd) {
    int N = c.path.size();       
    int len = N-1;               
    int m = 1 + (int)rnd.Rannyu(0, len/2);   		   // Block length
    int pos1 = 1 + (int)rnd.Rannyu(0, len - 2*m + 1);  // Start index of the first block
    int min_pos2 = pos1 + m;  						   // Second block must start after the first                    
    int max_pos2 = N - m;     						   // Last valid start index                      	
    int pos2 = min_pos2 + (int)rnd.Rannyu(0, max_pos2 - min_pos2 + 1);	// Start index of the second block

    for (int t=0; t<m; t++) {
		// Swap the two blocks element by element
        swap(c.path[pos1 + t], c.path[pos2 + t]);
    }
}

// ==================== CROSSOVER ====================

void crossover(Chromosome& child, Chromosome& parent2, int pos, int len) {
    int N = child.path.size();

    // Extract segment from child
    vector<int> block(len);
    for (int i=0; i<len; i++) block[i] = child.path[pos+i];

    // Rewrite the segment following order in parent2
    int count = 0;
    for (int i=1; i<N and count<len; i++) {
        int g = parent2.path[i];
        
        for (int j=0; j<len; j++) {
            if (block[j] == g) {
                child.path[pos + count] = g;
                count++;
                break;
            }
        }
    }
}

// ==================== EVOLUTION & METRICS ====================

void evolve_one_generation(Population& pop, vector<City>& cities, Random& rnd,
                           double psel, double pcross, double pm_swap, double pm_inv, double pm_shift, double pm_mperm) {
							   
    int M = pop.individuals.size();
    int N = pop.individuals[0].path.size();

    vector<Chromosome> next;
    next.reserve(M);
	
	// For each new individual
    for (int i=0; i<M; i++) {
		// Select two parents
        int i1, i2;
        select_two_parents(pop, rnd, psel, i1, i2);
        Chromosome& p1 = pop.individuals[i1];
        Chromosome& p2 = pop.individuals[i2];

		// Copy parent1 into child
        Chromosome child = p1;
		
		// Apply crossover with parent2, with probability pcross
        if (rnd.Rannyu()<pcross and N>2) {
            int pos = 1 + (int)rnd.Rannyu(0, N-1);    		// Segment start in [1..N-1]
            int len = 1 + (int)rnd.Rannyu(0, N-pos);  		// Segment length in [1..N-pos]
            crossover(child, p2, pos, len);
        }

		// Apply each mutation with given probability
        if (rnd.Rannyu() < pm_swap and N > 2) mutate_swap (child, rnd);
        if (rnd.Rannyu() < pm_inv and N > 2) mutate_inversion (child, rnd);
        if (rnd.Rannyu() < pm_shift and N > 2) mutate_shift(child, rnd);
        if (rnd.Rannyu() < pm_mperm and N > 4) mutate_mpermut (child, rnd);

		// Validate and re-evaluate fitness
        child.check();
        evaluate_fitness(cities, child);
		// Add child to new population
        next.push_back(child);
    }
	// Replace old population and sort
    pop.individuals.swap(next);
    sort_population(pop);
}

double best_length(Population& pop) {
	// Best individual at the head
    return pop.individuals[0].length;
}

double avg_best_half_length(Population& pop) {
    int M = pop.individuals.size();
    int half = max(1, M / 2);
    double s = 0.0;
    for (int i=0; i<half; i++) 
		s += pop.individuals[i].length;
    return s / half;
}

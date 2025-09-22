// EXERCISE 04_1 / 04_2

#include "system.h"

using namespace std;

int main (int argc, char *argv[]) {

	System SYS;
	SYS.initialize();							 // Read input parameters from 'input.dat'
	SYS.initialize_properties();				 // Read list of requested observables from 'properties.dat'
	SYS.block_reset(0);						 	 // Reset block averages
	
	int nconf = 1;
	
// Compute averages in each block
	for(int i=0; i < SYS.get_nbl(); i++){ 	 	 // Loop over blocks
	
		for(int j=0; j < SYS.get_nsteps(); j++){ // Loop over steps in a block
			SYS.step();							 // Perform a Monte Carlo move (Metropolis or Gibbs)
			SYS.measure();						 // Measure requested observables
	  
// Commented to avoid "filesystem full"! 
//      	if(j%50 == 0){
//        		SYS.write_XYZ(nconf); 			 // Write actual configuration in XYZ format
//        		nconf++;
//      	}
		}
		
		SYS.averages(i+1);						 // Compute and store block averages
		SYS.block_reset(i+1);					 // Reset accumulators for the next block
	}
	
	SYS.finalize();								 // Save final configuration

	return 0;
}

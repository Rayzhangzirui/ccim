// #define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
// #define CATCH_CONFIG_RUNNER
// #define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include "helper.h"
#include "tryvn.h"


//test connected components
TEST_CASE("write_surf"){
	GridData grid;
	int GRIDNUM = 50;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
	double radius = 0.5;
	// create initial surface
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	for(int opt = 11; opt<=16; opt++){
		init_surf(S, radius, grid, opt);
		write_field("initialsurface"+to_string(opt)+".txt", S, grid);	
	}	
}


#include <iostream>
#include <Solver.h>
#include <Thrust2DIncompressibleSolver.h>

#include <Grid.h>
#include <Fields.h>
#include <Params.h>
#include <Boundaries.h>

#include <read.h>


int main(){

	Grid   		grid;
	Fields 		fields;
	Params 		params;
	Boundaries  boundaries;

	read_inputs("simdata.yaml", grid, params, boundaries);
	std::cout << "OK" << std::endl;
	Thrust2DIncompressibleSolver solver(grid, fields, params, boundaries);

	std::cout << "initialising .. " << std::endl;
	solver.initialise();

	std::cout << "initialised" << std::endl;
	
	int i = 0;

	while(!solver.finished()){
		std::cout << i << std::endl;
		solver.take_step();
		i++;
	}

	solver.write_results();
	solver.close();

	return 0;
}

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

	std::cout << "Reading input data ... " << std::endl;
	read_inputs("../simdata.yaml", grid, params, boundaries);
	std::cout << "Initialising solver ... " << std::endl;
	Thrust2DIncompressibleSolver solver(grid, fields, params, boundaries);
	solver.initialise();
	
	std::cout << "Time steps: " << std::endl;
	int i = 0;
	while(!solver.finished()){
		std::cout << i << std::endl;
		solver.take_step();
		i++;
	}

	std::cout << "Writing results ..." << std::endl;
	solver.write_results();
	solver.close();

	return 0;
}

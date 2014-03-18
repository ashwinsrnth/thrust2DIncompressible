#include <iostream>
#include <Solver.h>

int main(){

	Grid   		grid;
	Fields 		fields;
	Params 		params;
	Boundaries  boundaries;

	Solver solver(grid, fields, params, boundaries);

	solver.initialise();

	while(!solver.finished()){
		solver.take_step();
	}

	solver.write_results();
	solver.close();

	return 0;
}
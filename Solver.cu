#include <Solver.h>


Solver::Solver(Grid& _g, Fields& _f, Params& _p, Boundaries& _b):
	grid(_g), params(_p), fields(_f), boundaries(_b){}

void Solver::initialise(){

}

void Solver::take_step(){

}

bool Solver::finished(){

	return true;
}

void Solver::write_results(){

}

void Solver::close(){

}



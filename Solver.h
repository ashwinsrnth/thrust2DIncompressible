#pragma once

#include <Grid.h>
#include <Fields.h>
#include <Params.h>
#include <Boundaries.h>

class Solver{
public:
	Solver(Grid&, Fields&, Params&, Boundaries&);
	void initialise();
	bool finished();
	void take_step();
	void write_results();
	void close();

private:
	Grid& grid;
	Fields& fields;
	Params& params;
	Boundaries& boundaries;
};
#pragma once

#include <Grid.h>
#include <Fields.h>
#include <Params.h>
#include <Boundaries.h>

class Solver{
/**
	Solver is a virtual base class for CFD solvers.
*/

public:
	virtual void initialise()=0;
	virtual bool finished()=0;
	virtual void take_step()=0;
	virtual void write_results()=0;
	virtual void close()=0;
};
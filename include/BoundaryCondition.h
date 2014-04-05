#pragma once
#include <types.h>


enum bc_t{
	DIRICHLET,
	NEUMANN,
};

class BoundaryCondition{
public:
	
	bc_t bc_type;
	Real bc_value;
};
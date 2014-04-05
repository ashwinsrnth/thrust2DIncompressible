#pragma once

#include <types.h>

class Fields{
public:
	// velocities, pressure and right-hand sides:
	Vector 	u, 
			v, 
			p, 
			r1x, 
			r1y, 
			r2;
};
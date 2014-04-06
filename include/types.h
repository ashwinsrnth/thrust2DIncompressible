#pragma once
#include <cusp/array1d.h>
#include <cusp/coo_matrix.h>

// Real is double or float:
typedef double Real;

// MemoryType is cusp::host_memory or cusp::device_memory:
typedef cusp::device_memory MemoryType;

// host and device vectors:
typedef cusp::array1d<Real, cusp::host_memory> 		HostVector;
typedef cusp::array1d<Real, cusp::device_memory>	DeviceVector;

// host and device int vectors:
typedef cusp::array1d<int,  cusp::host_memory> 		HostIntVector;
typedef cusp::array1d<int,  cusp::device_memory> 	DeviceIntVector;

// Vector types:
typedef cusp::array1d<Real, MemoryType>			    Vector;
typedef cusp::array1d<int,  MemoryType>			    IntVector;

// COO matrix and preconditioner:
typedef cusp::coo_matrix<int, Real, MemoryType> 	COOMatrix;


#pragma once
#include <solvers/Solver.h>
#include <types.h>



class Thrust2DIncompressibleSolver : public Solver{

public:
	Thrust2DIncompressibleSolver(Grid&, Fields&, Params&, Boundaries&);
	void initialise();
	bool finished();
	void take_step();
	void write_results();
	void close();

private:
	void initialise_arrays();
	void make_labels();	
	void make_poisson_matrix();
	void compute_rhs1();
	void compute_intermediate_velocity();
	void update_ghosts();
	void compute_rhs2();
	void solve_poisson_system();
	void update_velocity();

	// storage for labels:
	IntVector u_labels;
	IntVector v_labels;
	IntVector p_labels;
	
	// a quirk is that we also need a host copy for 
	// p_labels:
	HostIntVector   p_labels_h;

	// I, J, V and COO matrix for Poisson system:
	HostIntVector I;
	HostIntVector J;
	HostVector    V;
	COOMatrix pMat;
	// number of non-zero entries in Poisson matrix:
	int num_entries;

	// number of solver steps:
	int num_steps;

	// simulation data:
	Grid& grid;
	Fields& fields;
	Params& params;
	Boundaries& boundaries;
};

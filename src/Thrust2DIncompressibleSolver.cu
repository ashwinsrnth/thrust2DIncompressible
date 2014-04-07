	#include <types.h>

	#include <Thrust2DIncompressibleSolver.h>
	#include <thrust_iterators.h>
	#include <functors.h>
	#include <write.h>

	#include <thrust/gather.h>
	#include <thrust/tuple.h>
	#include <thrust/for_each.h>
	#include <thrust/fill.h>

	#include <cusp/krylov/bicgstab.h>
	#include <cusp/krylov/cg.h>

	Thrust2DIncompressibleSolver::Thrust2DIncompressibleSolver(Grid& _g, Fields& _f, Params& _p, Boundaries& _b):
		grid(_g),
		fields(_f),
		params(_p),
		boundaries(_b){}

	void Thrust2DIncompressibleSolver::initialise(){
		initialise_arrays();
		make_labels();
		make_poisson_matrix();
	}

	void Thrust2DIncompressibleSolver::take_step(){
		update_ghosts();
		compute_rhs1();
		compute_intermediate_velocity();
		update_ghosts();
		compute_rhs2();
		solve_poisson_system(); 
		update_velocity();
		num_steps += 1;
	}

	bool Thrust2DIncompressibleSolver::finished(){
		return (num_steps >= params.nsteps);
	}

	// TODO
	void Thrust2DIncompressibleSolver::write_results(){

		HostVector u_h = fields.u;
		HostVector v_h = fields.v;
		HostVector p_h = fields.p;

		write_vector<Real>(u_h, "u.txt");
		write_vector<Real>(v_h, "v.txt");
		write_vector<Real>(p_h, "p.txt");

	}

	// TODO
	void Thrust2DIncompressibleSolver::close(){

	}



	//======================================================================//
	//								INITIALISE 								//
	//======================================================================//

	void Thrust2DIncompressibleSolver::initialise_arrays(){

		int N = grid.N_x*grid.N_y;

		fields.u.resize(N, 0); 		
		fields.v.resize(N, 0);			
		fields.r1x.resize(N, 0);		
		fields.r1y.resize(N, 0);		
		fields.r2.resize(N, 0);		
		fields.p.resize(N, 0);			
		u_labels.resize(N, 0);			
		v_labels.resize(N, 0);			
		p_labels.resize(N, 0);			
		p_labels_h.resize(N, 0);	

		// For solving Poisson system:	
		num_entries = 	(grid.N_x - 4)*(grid.N_y - 4)*5 + 	// inner
						4*3 +								// corners
						2*(grid.N_x-4)*4 +					// edges 
						2*(grid.N_y-4)*4 +					// edges
						2*grid.N_x + 2*(grid.N_y-2);		// outer

		I.resize(5*N, 0);	
		J.resize(5*N, 0); 	
		V.resize(5*N, 0);
	}

	void Thrust2DIncompressibleSolver::make_labels(){
		thrust::for_each_n(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						u_labels.begin(),
						v_labels.begin(),
						p_labels.begin(),
						thrust::make_counting_iterator(0))),
				grid.N_x*grid.N_y,
				make_labels_functor(grid));
	}


	void Thrust2DIncompressibleSolver::make_poisson_matrix(){
		std::cout << p_labels.size() << std::endl;
		std::cout << p_labels_h.size() << std::endl;
		
		p_labels_h = p_labels;

		thrust::for_each_n(	thrust::make_zip_iterator(
								thrust::make_tuple(
									thrust::make_counting_iterator(0),
									thrust::make_zip_iterator(
										thrust::make_tuple(
											p_labels_h.begin(),
											p_labels_h.begin()-1,
											p_labels_h.begin()+1,
											p_labels_h.begin()-grid.N_x,
											p_labels_h.begin()+grid.N_x)))),
							grid.N_x*grid.N_y,
							make_poisson_functor(I, J, V, grid));


		thrust::remove_if(	I.begin(),
								I.end(),
								V.begin(),
								is_zero());
		thrust::remove_if(	J.begin(),
								J.end(),
								V.begin(),
								is_zero());
		thrust::remove_if(	V.begin(),
								V.end(),
								V.begin(),
								is_zero());

		I.resize(num_entries); 
		J.resize(num_entries); 
		V.resize(num_entries);

		pMat.resize(grid.N_x*grid.N_y, grid.N_x*grid.N_y, num_entries);
		pMat.row_indices = I;
		pMat.column_indices = J;
		pMat.values = V; 
		pMat.sort_by_row_and_column();

		/* preconditioner: */
        preconditioner = new cusp::precond::aggregation::smoothed_aggregation<int, Real, MemoryType>(pMat);
 
	}
	

	//======================================================================//
	//							   SOLVER STEP 								//
	//======================================================================//


	void Thrust2DIncompressibleSolver::compute_rhs1(){

		thrust::for_each_n(
				thrust::make_zip_iterator(
					thrust::make_tuple(
					fields.r1x.begin(),
					fields.r1y.begin(),
					u_labels.begin(),
					v_labels.begin(),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							fields.u.begin(),
							fields.u.begin()-1,
							fields.u.begin()+1,
							fields.u.begin()-grid.N_x,
							fields.u.begin()+grid.N_x)),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							fields.v.begin()-grid.N_x,
							fields.v.begin()-grid.N_x+1,
							fields.v.begin(),
							fields.v.begin()+1)),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							fields.v.begin(),
							fields.v.begin()-1,
							fields.v.begin()+1,
							fields.v.begin()-grid.N_x,
							fields.v.begin()+grid.N_x)),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							fields.u.begin()-1,
							fields.u.begin(),
							fields.u.begin()+grid.N_x-1,
							fields.u.begin()+grid.N_x)),
					thrust::make_counting_iterator(0))),
				grid.N_x*grid.N_y,
				compute_rhs1_functor(grid, params));
	}

	void Thrust2DIncompressibleSolver::compute_intermediate_velocity(){

		// TODO: Use cusp::blas functions instead of custom saxpy functor:

		thrust::for_each_n(thrust::make_zip_iterator(
							thrust::make_tuple(
								fields.u.begin(),
								fields.u.begin(),
								fields.r1x.begin())),
							grid.N_x*grid.N_y,
							saxpy_functor(params.dt));

		thrust::for_each_n(thrust::make_zip_iterator(
							thrust::make_tuple(
								fields.v.begin(),
								fields.v.begin(),
								fields.r1y.begin())),
							grid.N_x*grid.N_y,
							saxpy_functor(params.dt));

	}

	void Thrust2DIncompressibleSolver::compute_rhs2(){
		thrust::for_each_n(thrust::make_zip_iterator(
							thrust::make_tuple(
								fields.r2.begin(),
								fields.u.begin(),
								fields.u.begin()-1,
								fields.v.begin(),
								fields.v.begin()-grid.N_x,
								p_labels.begin(),
								thrust::make_counting_iterator(0))),
							grid.N_x*grid.N_y,
							compute_rhs2_functor(grid, params));

	}


	void Thrust2DIncompressibleSolver::solve_poisson_system(){
	    // cusp::verbose_monitor<double> monitor(fields.r2, 500, 1e-3);
	    cusp::default_monitor<double> monitor(fields.r2, 500, 1e-3);
	    // cusp::identity_operator<double, cusp::device_memory> M(pMat.num_rows, pMat.num_rows);
	    // cusp::krylov::bicgstab(pMat, fields.p, fields.r2, monitor);
	    // cusp::krylov::cg(pMat, fields.p, fields.r2, monitor, *preconditioner);
	    cusp::krylov::bicgstab(pMat, fields.p, fields.r2, monitor, *preconditioner);
	}


	void Thrust2DIncompressibleSolver::update_velocity(){
		
		thrust::for_each_n(thrust::make_zip_iterator(
							thrust::make_tuple(
								fields.u.begin(),
								fields.v.begin(),
								thrust::make_zip_iterator(
									thrust::make_tuple(
										fields.p.begin(),
										fields.p.begin()+1,
										fields.p.begin()+grid.N_x)),
								u_labels.begin(),
								v_labels.begin())),
							grid.N_x*grid.N_y,
							update_velocity_functor(grid, params));
	}



	void Thrust2DIncompressibleSolver::update_ghosts(){

		thrust::for_each_n(thrust::make_zip_iterator(
							thrust::make_tuple(
								fields.u.begin(),
								fields.u.begin()-1,
								fields.u.begin()+1,
								fields.u.begin()-grid.N_x,
								fields.u.begin()+grid.N_x,
								u_labels.begin(),
								thrust::make_counting_iterator(0))),
						grid.N_x*grid.N_y,
						update_ghosts_functor(grid, boundaries));

		thrust::for_each_n(thrust::make_zip_iterator(
								thrust::make_tuple(
								fields.v.begin(),
								fields.v.begin()-1,
								fields.v.begin()+1,
								fields.v.begin()-grid.N_x,
								fields.v.begin()+grid.N_x,
								v_labels.begin(),
								thrust::make_counting_iterator(0))),
						grid.N_x*grid.N_y,
						update_ghosts_functor(grid, boundaries));

	}

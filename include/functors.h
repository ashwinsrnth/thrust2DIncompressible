#include <types.h>

/*
	Functors defined here:

	Compute functors:
		make_labels_functor
		calculate_intermediate_velocity_functor
		make_poisson_functor

	Conditional functors:
		is_ghost
*/

class is_ghost{

public:
	__host__ __device__
	bool operator () (int x){
		if (x == 2)
			return true;
	
	return false;
	}
};


class is_zero{
public:
	bool operator () (Real x){
		if (x == 0){
			return true;
		}
		else
			return false;
	}
};


class make_labels_functor{

	/**
		Attach labels for velocity and 
		ressure fields.
		This is one way to tell the solver which points
		are compute, boundary, ghost or padding. For
		velocity fields (u and v):

		0 - compute, 1 - boundary, 2 - ghost, 3 - padding

		For pressure field, the only labels are

		0 - inner, 1 - outer
	
		Parameters to functor:

		thrust::Tuple T

		T is dereferenced to get:
		<u_label, v_label, p_label, thrust::counting_iterator>

		The counting iterator gives the index of the point.
	*/

public:
	make_labels_functor(Grid _grid):grid(_grid){}

	
	template <typename Tuple>
	__host__ __device__
	void operator () (Tuple T){

		ix = thrust::get<3>(T)%grid.N_x;
		iy = thrust::get<3>(T)/grid.N_x;

		// u_labels

			// boundary:
			if (ix == 0 or ix == grid.N_x - 2){
				thrust::get<0>(T) = 1;
			}

			// ghost:
			else if ((iy == 0 or iy == grid.N_y - 1) and (ix > 0 and ix < grid.N_x - 1)){
				thrust::get<0>(T) = 2;
			}

			// pad:

			else if (ix == grid.N_x-1){
				thrust::get<0>(T) = 3;
			}

		// v_labels

			// boundary:
			if (iy == 0 or iy == grid.N_y - 2){
				thrust::get<1>(T) = 1;
			}

			// ghost:
			else if ((ix == 0 or ix == grid.N_x - 1) and (iy > 0 and iy < grid.N_x - 1)){
				thrust::get<1>(T) = 2;
			}

			//

			else if (iy == grid.N_y - 1){
				thrust::get<1>(T) = 3;
			}
	
		// p_labels

			if (ix == 0 or iy == 0 or ix == grid.N_x - 1 or iy == grid.N_y - 1)
				thrust::get<2>(T) = 1;
	}

private:
	Grid grid;
	int ix, iy;

};


class saxpy_functor{
public:
	saxpy_functor(Real _x):x(_x){}

	template <typename Tuple>
	__host__ __device__
	void operator () (Tuple T){
		thrust::get<0>(T) = thrust::get<1>(T) + x*thrust::get<2>(T);
	}

private:
	Real x;
};



class compute_rhs1_functor{

	/**

		Parameters to functor:

		thrust::tuple T

		T is dereferenced to get:

		<u_star, v_star, u_labels, v_labels, T1, T2, T3, T4, thrust::counting_iterator>
		
		u_star and v_star are the intermediate velocities.

		T1 -> thrust::tuple<u, u-1, u+1, u-N_x, u+N_x> // stencil for
		T2 -> thrust::tuple<v-N_x, v-N_x+1, v, v+1>	   // updating u
		T3 -> thrust::tuple<v, v-1, v+1, v-N_x, v+N_x> // stencil for
		T4 -> thrust::tuple<u-1, u, u+N_x-1, u+N_x>    // updating v

		The counting iterator gives the index of the point.
	*/

public:

	compute_rhs1_functor(const Grid _grid, const Params _params):
										grid(_grid),
										params(_params){
		dx = grid.L_x/grid.N_x;
		dy = grid.L_y/grid.N_y;
	}


	template <typename Tuple>
	__host__ __device__
	void operator () (Tuple T){

		// x- and y- indices:
		ix = thrust::get<8>(T)%grid.N_x;
		iy = thrust::get<8>(T)/grid.N_x;

		// compute u_star and v_star:

		if (thrust::get<2>(T) == 0){

			// Aliases:
			u 	= thrust::get<0>(thrust::get<4>(T));
			ul 	= thrust::get<1>(thrust::get<4>(T));
			ur  = thrust::get<2>(thrust::get<4>(T));
			ud  = thrust::get<3>(thrust::get<4>(T));
			uu  = thrust::get<4>(thrust::get<4>(T));
			vll = thrust::get<0>(thrust::get<5>(T));
			vlr = thrust::get<1>(thrust::get<5>(T));
			vul = thrust::get<2>(thrust::get<5>(T));
			vur = thrust::get<3>(thrust::get<5>(T));

			// compute r1x
			thrust::get<0>(T) 	=
									((u + ul)*(u + ul) - (u + ur)*(u + ur))/(4.*dx) +
									((u + ud)*(vll + vlr) - (u + uu)*(vul + vur))/(4.*dy) +
									(1./params.Re)*
									((ul - 2*u + ur)/(dx*dx) + (ud - 2*u + uu)/(dy*dy));
		}

		if (thrust::get<3>(T) == 0){

			// Aliases
			v 	= thrust::get<0>(thrust::get<6>(T));
			vl 	= thrust::get<1>(thrust::get<6>(T));
			vr  = thrust::get<2>(thrust::get<6>(T));
			vd  = thrust::get<3>(thrust::get<6>(T));
			vu  = thrust::get<4>(thrust::get<6>(T));
			ull = thrust::get<0>(thrust::get<7>(T));
			ulr = thrust::get<1>(thrust::get<7>(T));
			uul = thrust::get<2>(thrust::get<7>(T));
			uur = thrust::get<3>(thrust::get<7>(T));

			// compute r1y
			thrust::get<1>(T) = 
									((v + vd)*(v + vd) - (v + vu)*(v + vu))/(4.*dy) +
									((v + vl)*(ull + uul) - (v + vr)*(ulr + uur))/(4.*dx) +
									(1./params.Re)*
									((vl - 2*v + vr)/(dx*dx) + (vd - 2*v + vu)/(dy*dy));
		}

	}


private:

	int ix, iy;
	Real dx, dy;

	// left, right, down, up, upper left, upper right, lower left, lower right
	Real u, ul, ur, ud, uu, vll, vlr, vul, vur;
	Real v, vl, vr, vd, vu, ull, ulr, uul, uur;
	const Grid grid;
	const Params params;
};


class make_poisson_functor{

/**
	Construct the 2-D sparse pressure Poisson matrix (host only).
	
	Parameters to constructor:
	(I, J, V, params)
	I, J and V are HostVectors for the row, column
	and values of the COO sparse matrix.

	Parameters to functor:
	thrust::tuple<thrust::counting_iterator(0), 
				  thrust::tuple<p, p_left, p_right, p_bottom, p_top> >,
	where p are the pressure labels.
*/

public:
	make_poisson_functor(HostIntVector& _I,
						 HostIntVector& _J,
						 HostVector&	_V,
						 Grid&			_grid):
	
						I(_I),
						J(_J),
						V(_V),
						grid(_grid){
	
		dx = grid.L_x/grid.N_x;
		dy = grid.L_y/grid.N_y;
	}

	template <typename Tuple>
	__host__
	void operator () (Tuple T){

		idx  = thrust::get<0>(T);
		i = 5*idx;

		if ((thrust::get<0>(thrust::get<1>(T)) == 1)){
			I[i] = idx; J[i] = idx - grid.N_x; 	V[i] = 0; i += 1;
			I[i] = idx; J[i] = idx - 1;	 		V[i] = 0; i += 1;
			I[i] = idx; J[i] = idx;		 		V[i] = 1; i += 1;
			I[i] = idx; J[i] = idx + 1;	 		V[i] = 0; i += 1;
			I[i] = idx; J[i] = idx + grid.N_x; 	V[i] = 0;
		}

		// if inner:
		else {

			diag = 0; left = 0; right = 0; bottom = 0; top = 0;

			if (thrust::get<1>(thrust::get<1>(T)) == 0){
				left  	= 1/(dx*dx);
				diag   -= 1/(dx*dx);
			}

			if (thrust::get<2>(thrust::get<1>(T)) == 0){
				right	= 1/(dx*dx);
				diag   -= 1/(dx*dx);
			}

			if (thrust::get<3>(thrust::get<1>(T)) == 0){
				bottom  = 1/(dy*dy);
				diag   -= 1/(dy*dy);
			}

			if (thrust::get<4>(thrust::get<1>(T)) == 0){
				top 	= 1/(dy*dy);
				diag   -= 1/(dy*dy);
			}

			I[i] = idx; J[i] = idx - grid.N_x; 	V[i] = bottom;	i += 1;
			I[i] = idx; J[i] = idx - 1;	 		V[i] = left; 	i += 1;
			I[i] = idx; J[i] = idx;		 		V[i] = diag;	i += 1;
			I[i] = idx; J[i] = idx + 1;	 		V[i] = right;   i += 1;
			I[i] = idx; J[i] = idx + grid.N_x; 	V[i] = top;

			if (idx == grid.N_x*(grid.N_y/2) + grid.N_x/2){
				// retrace and correct:
				V[i-2] += 1;
			}
		}	



	}

	HostIntVector& 	I;	
	HostIntVector& 	J;
	HostVector&   	V;

private:
	Grid& grid;
	Real dx, dy;
	int i;
	int idx;
	int diag, left, right, bottom, top;

};


class compute_rhs2_functor {

public:

	compute_rhs2_functor(Grid _grid, Params _params):grid(_grid), params(_params){
		dx = grid.L_x/grid.N_x;
		dy = grid.L_y/grid.N_y;
	}

	template <typename Tuple>
	__host__ __device__

	void operator () (Tuple T){
		if (thrust::get<5>(T) == 0){
			if (thrust::get<6>(T) == grid.N_x*grid.N_x/2 + grid.N_x/2)
				thrust::get<0>(T) = 0;
			else{
			thrust::get<0>(T) = ((thrust::get<1>(T) - thrust::get<2>(T))/dx + 
								 (thrust::get<3>(T) - thrust::get<4>(T))/dy)/params.dt;
			}
		}
	}

private:
	Real dx, dy;
	Grid grid;
	Params params;
};


class update_velocity_functor {

public:

	update_velocity_functor(Grid _grid, Params _params):grid(_grid),params(_params){
		dx = grid.L_x/grid.N_x;
		dy = grid.L_y/grid.N_y;
	}

	template <typename Tuple>
	__host__ __device__
	void operator () (Tuple T){
		if (thrust::get<3>(T) == 0){
			thrust::get<0>(T) -=
								params.dt*(thrust::get<1>(thrust::get<2>(T)) -
										thrust::get<0>(thrust::get<2>(T)))/dx;
		}

		if (thrust::get<4>(T) == 0){
			thrust::get<1>(T) -=
								params.dt*(thrust::get<2>(thrust::get<2>(T)) -
										thrust::get<0>(thrust::get<2>(T)))/dy;
		}
	}

private:
	Grid grid;
	Params params;
	Real dx, dy;
};




class update_ghosts_functor{

/**
	Tuple:
	<V, V-1, V+1, V-N_X, V+N_x, counting_iterator>

*/

public:
	update_ghosts_functor(Grid _grid, Boundaries _boundaries):grid(_grid), boundaries(_boundaries){}

	template <typename Tuple>
	__host__ __device__
	void operator () (Tuple T){

		if (thrust::get<5>(T) == 2){

			ix = thrust::get<6>(T)%grid.N_x;
			iy = thrust::get<6>(T)/grid.N_x;

			if (ix == 0){
				// left
				if (boundaries.left.bc_type == DIRICHLET){
					thrust::get<0>(T) = 2*boundaries.left.bc_value - thrust::get<2>(T);
				}
			}

			else if (ix == grid.N_x-1){
				// right
				if (boundaries.right.bc_type == DIRICHLET){
					thrust::get<0>(T) = 2*boundaries.right.bc_value - thrust::get<1>(T);
				}
			}

			else if (iy == 0){
				// bottom
				if (boundaries.bottom.bc_type == DIRICHLET){
					thrust::get<0>(T) = 2*boundaries.bottom.bc_value - thrust::get<4>(T);
				}
			}

			else if (iy == grid.N_y-1){
				// top
				if (boundaries.top.bc_type == DIRICHLET){
					thrust::get<0>(T) = 2*boundaries.top.bc_value - thrust::get<3>(T);
				}
			}
		}
	}

private:
	Grid grid;
	Boundaries boundaries;
	int ix, iy;

};






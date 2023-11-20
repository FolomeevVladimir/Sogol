#pragma once
#include "mpi.h"
#include <cstdio>
#include <iostream>
#include "Vectorxd.h"

namespace sogol {

	//template<int d=1>
	class HPC {
		public:  // temporaly public
		static const int d = 1;
		MPI::Cartcomm cartcomm;
		int mpi_size, mpi_rank, rank, maxdims, source, up, down;
		int dims[d];
		int coords[d], coords_up[d], coords_down[d];
		bool periods[d];
		float tmp;
		MPI::Status state;
	public:
		HPC(sogol::Vectorxd<d> per = {1}) {
			MPI::Init();
			mpi_size = MPI::COMM_WORLD.Get_size();
			mpi_rank = MPI::COMM_WORLD.Get_rank();
			for (int i = 0; i < d; i++) {
				dims[i] = mpi_size;
				periods[i] = per[i];
			
			}
			
			cartcomm = MPI::COMM_WORLD.Create_cart(d, dims, periods, false);
			cartcomm.Get_coords(mpi_rank, d, coords);
			for (int i = 0; i < d; i++) {
				coords_up[i] = coords_down[i] = coords[i];
			}
			coords_up[0]++;
			coords_down[0]--;
			up = cartcomm.Get_cart_rank(coords_up);
			down = cartcomm.Get_cart_rank(coords_down);

				
		}
		void swap(float* f,float* nextf) {
			tmp = *nextf;
			cartcomm.Send(f, 1, MPI::FLOAT, up, 1);
			cartcomm.Recv(nextf, 1, MPI::FLOAT, down, 1, state);

			cartcomm.Send(&tmp, 1, MPI::FLOAT, down, 2);
			cartcomm.Recv(f, 1, MPI::FLOAT, up, 2, state);
			
			
		}
		void exit() {
			MPI::Finalize();
		}

	};











}
#pragma once
#include <iostream>
#include "Vectorxd.h"
#include <cstdlib>
#include <vector>
#include "Cell.h"
#include "Lattice.h"
#include "Stencil.h"
#include "Physics.h"
#include <fstream>
#include <chrono>
#include "BC.h"
#include "HPC.h"

int main() {



	//auto x= sogol::Vectorxd<3, int>();
	auto x = sogol::Vectorxd<2, double>();
	using Vector2d = sogol::Vectorxd<2, double>;
	Vector2d y = Vector2d();
	using d2q9 = sogol::D2Q9;
	const unsigned d = d2q9::d, q = d2q9::q;
	sogol::Vectorxd<d, unsigned> L = { 5,5 };
	sogol::HPC hpc;
	sogol::Lattice<d2q9> grid(L,&hpc);

	
	


	sogol::BGK<d2q9> bgk;
	sogol::Cell<2, 9> cell;



	grid.init(bgk, { 0.2,0. }, 1);
	sogol::Settype<2, 9> set(L);
	sogol::F<2, 9> f;


	std::cout<<cell<< std::endl;
	
	auto start = std::chrono::system_clock::now();
	for (unsigned i = 0; i < 1000; i++) {
		//grid.collide_stream(bgk);
		grid.collide(bgk);
		grid.stream();
		grid.applyBC();
		
		
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> sec = end - start;
	std::cout << "time sec=" << sec.count() << std::endl;

	sogol::Save<2, 9> save;
	save.init(L, hpc.mpi_rank);
	grid.iterate(save);




	double t1 = -1, t2 = hpc.mpi_rank;
	hpc.swap(&t1, &t2);
	std::cout << hpc.mpi_rank << "::" << t1 << "	" << t2 << "\n";

	
	if (hpc.mpi_rank == 0) {
		grid.cuInit();
		grid.cuRun();
		}
	hpc.exit();

	
	


}

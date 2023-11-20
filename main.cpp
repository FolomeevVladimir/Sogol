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
	auto x = sogol::Vectorxd<2, float>();
	using Vector2d = sogol::Vectorxd<2, float>;
	Vector2d y = Vector2d();
	using d2q9 = sogol::DQ<2,9>;
	const int d = d2q9::d, q = d2q9::q;
	sogol::Vectorxd<d, int> L = { 10,10 };
	sogol::HPC hpc;
	

	
	


	using bgk=sogol::BGK<d2q9>;
	sogol::Cell<2, 9> cell;

	sogol::Lattice<d2q9, sogol::BGK<d2q9> > grid(L, &hpc);

	grid.init( { 0.1,0. }, 1);
	sogol::Settype<2, 9> set(L);
	sogol::F<2, 9> f;


	std::cout<<cell<< std::endl;
	
	auto start = std::chrono::system_clock::now();
	grid.cuInit();
	for (int i = 0; i < 10; i++) {
		
		
		grid.cuRun();
		grid.stream();
		
		grid.applyBC();
		std::cout << "iteration=" << i << std::endl;
		
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<float> sec = end - start;
	std::cout << "time sec=" << sec.count() << std::endl;
	sogol::Save<2, 9> save;
	save.init(L, hpc.mpi_rank);
	grid.iterate(save);

	//std::cout << grid;
	



	
	
	
	//hpc.exit();

	
	


}

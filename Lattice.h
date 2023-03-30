#pragma once
#include <cmath>
#include <iostream>
#include <assert.h> 
#include "Vectorxd.h"
#include "Cell.h"
#include "Functors.h"
#include "Stencil.h"
#include "Physics.h"
#include <omp.h>
#include "BC.h"
#include <functional>
#include "HPC.h"
#include "Engine.cuh"

namespace sogol {

template<typename dq>
class Lattice {
public:
	static constexpr unsigned d = dq::d, q = dq::q;
	static constexpr unsigned half = (q - 1) / 2;
	Cell<d,q> *c;
	Vectorxd<d, unsigned> mask;
	Vectorxd<d, unsigned> dims;
	Vectorxd<d, unsigned> ind;
	std::vector<BC<dq>*> bcs;
	int size;
	sogol::HPC *hpc;
	//obsolete
	inline void do_for(unsigned k, func<d,q> &f) {

		if (k == d) {
			f( (*this)(ind),ind);
			return;
		}
		for (unsigned i = 0; i < dims[k]; i++) {
			ind[k] = i;
			do_for(k + 1,f);
		}
	}
	
public:
	inline Cell<d,q>& operator[] (unsigned i) {
		assert(i < size);
		return c[i];
	}
	Lattice(Vectorxd<d, unsigned int> L,sogol::HPC* hp) {
		dims = L;
		size = dims.conv();
		c = new Cell<d,q>[size];
		mask = {};
		ind = {};
		hpc = hp;
		
		for (unsigned i = 0; i <d; i++) {
			unsigned res = 1;
			for (unsigned j = i + 1; j < d; j++) {
				res *= L[j];				
			}
			mask[i] = res;
		}
		std::cout << mask;
	}
	inline Cell<d, q>& operator() (Vectorxd<d, unsigned> p) {
		assert((p * mask) < size);
		return c[(p * mask)];
	}
	//obsolete
	void iterate(func<d,q> &f) {
		
		do_for(0,f);
	}
	void spatial_for(std::function<void(sogol::Cell<d, q> & c, sogol::Vectorxd<d, unsigned> ind)> f) {

		_for(0, f);
	}
	inline void _for(unsigned k, std::function<void(sogol::Cell<d, q> & c, sogol::Vectorxd<d, unsigned> ind)> f) {

		if (k == d) {
			f((*this)(ind), ind);
			return;
		}
		for (unsigned i = 0; i < dims[k]; i++) {
			ind[k] = i;
			_for(k + 1, f);
		}
	}
	void init(Dynamics<dq> &p,Vectorxd<d,double> u,double rho) {
		for (int i = 0; i < size; i++) {
			c[i].f = p.feq(u, rho);
		
		}
		auto per = new MPIexchange<dq>(this->mask, this->c, 3, { (int)dims[0],0},hpc);
		bcs.push_back(per);
		auto per2 = new MPIexchange<dq>(this->mask, this->c, 4, { -1*(int)dims[0],0},hpc);
		bcs.push_back(per2);
		std::cout << "period " << ((Periodic<dq>* )bcs[1])->period << std::endl;
		
		spatial_for([&](sogol::Cell<d, q>& c, sogol::Vectorxd<d, unsigned> ind) mutable{

			if (ind[1] == dims[1] - 1) { c.type = 2; }
			if (ind[1] == 0) { c.type = 2; }

			if (ind[0] == 0 && ind[1] != 0 && ind[1] != (dims[1] - 1)) { c.type = 3; bcs[0]->indexes.push_back(ind * mask); }
			if (ind[0] == dims[0] - 1 && ind[1] != 0 && ind[1] != (dims[1] - 1)) { c.type = 4; bcs[1]->indexes.push_back(ind * mask);}
			});
	

	}
	void collide(Dynamics<dq>& p) {
		Vectorxd<d, double> u;
		double rho;
		Vectorxd<q, double> feq;

#pragma omp parallel for
		for (int i = 0; i < size; i++) {

			

			if (c[i].type != 2) {

				p.collide(c[i]);

				c[i].swap();
			}
	
		}
	}
	void stream() {
	
		int next;
		
		auto ck = dq::c;
#pragma omp parallel for private(next)
		for (int i = 0; i < size; i++) {
			if (c[i].type == 1) {
				for (unsigned k = 1; k <= half; k++) {
					next = i + ck[k] * mask;
					if (next < 0) { continue; }
					if (next >= size) { continue; }

				

					assert(next >= 0);
					assert(next < size);
					std::swap(c[i].opposite(k), c[next][k]);


				}
			}
		}
		
	
	}

	void collide_stream(Dynamics<dq>& p) {
		Vectorxd<d, double> u;
		double rho;
		Vectorxd<q, double> feq;
		auto ck = dq::c;

		//#pragma omp parallel for
		for (int i = 0; i < size; i++) {

			if (c[i].type == 1) {

				p.collide(c[i]);

				//c[i].swap();
			}
			if (c[i].type == 3) {

				p.collide(c[i]);

				c[i].swap();
			}
			if (c[i].type == 4) {

				p.collide(c[i]);

				c[i].swap();
			}
			int next ;
			for (unsigned k = 1; k <= half; k++) {
				if (c[i].type !=1 ) { break; }
				next = i + ck[k] * mask;
				if (next < 0) { continue; }
				if (next >= size) { continue; }
				
				assert(next >= 0);
				assert(next < size);
				
				upwind_swap(next, i, k);
				}
			
		}
		

}
	inline void applyBC() {
		
		for (auto b : bcs) {
			b->process();
		}
	}
	inline void upwind_swap(int next, int i, unsigned k) {
			double ftmp = c[i][k];
			c[i][k] = c[i].opposite(k);
			c[i].opposite(k) = c[next][k];
			c[next][k] = ftmp;
		}
	void cuInit() {
		if (hpc->mpi_rank == 0){
			std::cout << "call initCu="<<(*hpc).mpi_rank<<"\n";
			
			auto cell=initCu<d,q>(c,size);
			delete[] c;
			c = cell;
		

		}

	}
	void cuRun() {
		if (hpc->mpi_rank == 0) {
			
			runCu(c, size);
			for (int i = 0; i < size; i++) {
				std::cout << c[i] << "	";
			}

		}

	}


	~Lattice() {
		delete[] c;
		
	}
	

};






}
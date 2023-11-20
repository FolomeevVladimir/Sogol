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

#include <concepts>



template <typename T>
concept Stencil =
	requires(T s) {
		{s.d}->std::same_as<int>;
		{s.q}->std::same_as<int>;
		{s.cs2}->std::same_as <float>;


};

namespace sogol {

template<typename dq,typename phys>
class Lattice {
public:
	static constexpr int d = dq::d, q = dq::q;
	static constexpr int half = (q - 1) / 2;
	Cell<d,q> *c;
	Vectorxd<d, int> mask;
	Vectorxd<d, int> dims;
	Vectorxd<d, int> ind;
	std::vector<BC<dq>*> bcs;
	int size;
	sogol::HPC *hpc;
	//obsolete
	inline void do_for(int k, func<d,q> &f) {

		if (k == d) {
			f( (*this)(ind),ind);
			return;
		}
		for (int i = 0; i < dims[k]; i++) {
			ind[k] = i;
			do_for(k + 1,f);
		}
	}
	
public:
	inline Cell<d,q>& operator[] (int i) {
		assert(i < size);
		return c[i];
	}
	Lattice(Vectorxd<d,  int> L,sogol::HPC* hp) {
		std::cout << "size=" << L;
		dims = L;
		size = dims.conv();
		std::cout << "size=" << size;
		c = new Cell<d,q>[size];
		mask = {};
		ind = {};
		hpc = hp;
		
		for (int i = 0; i <d; i++) {
			int res = 1;
			for (int j = i + 1; j < d; j++) {
				res *= L[j];				
			}
			mask[i] = res;
		}
		std::cout << mask;
	}
	inline Cell<d, q>& operator() (Vectorxd<d, int> p) {
		assert((p * mask) < size);
		return c[(p * mask)];
	}
	//obsolete
	void iterate(func<d,q> &f) {
		
		do_for(0,f);
	}
	void spatial_for(std::function<void(sogol::Cell<d, q> & c, sogol::Vectorxd<d, int> ind)> f) {

		_for(0, f);
	}
	inline void _for(int k, std::function<void(sogol::Cell<d, q> & c, sogol::Vectorxd<d, int> ind)> f) {

		if (k == d) {
			f((*this)(ind), ind);
			return;
		}
		for (int i = 0; i < dims[k]; i++) {
			ind[k] = i;
			_for(k + 1, f);
		}
	}
	void init(Vectorxd<d,float> u,float rho) {
		for (int i = 0; i < size; i++) {
			c[i].f =  phys::feq(u, rho);
		
		}
		auto per = new MPIexchange<dq>(this->mask, this->c, 3, { (int)dims[0],0 }, hpc);
		bcs.push_back(per);
		auto per2 = new MPIexchange<dq>(this->mask, this->c, 4, { -1 * (int)dims[0],0 }, hpc);
		bcs.push_back(per2);
		std::cout << "period " << ((Periodic<dq>* )bcs[1])->period << std::endl;
		
		spatial_for([&](sogol::Cell<d, q>& c, sogol::Vectorxd<d, int> ind) mutable{

			if (ind[1] == dims[1] - 1) { c.type = 2; }
			if (ind[1] == 0) { c.type = 2; }

			if (ind[0] == 0 && ind[1] != 0 && ind[1] != (dims[1] - 1)) { c.type = 3; bcs[0]->indexes.push_back(ind * mask); }
			if (ind[0] == dims[0] - 1 && ind[1] != 0 && ind[1] != (dims[1] - 1)) { c.type = 4; bcs[1]->indexes.push_back(ind * mask);}
			});
	

	}
	void collide() {
		Vectorxd<d, float> u;
		float rho;
		Vectorxd<q, float> feq;
		std::cout << "size=" << size;
#pragma omp parallel for
		for (int i = 0; i < size; i++) {

			

			if (c[i].type != 2) {

				phys::collide(c[i]);

				c[i].swap();
			}
	
		}
	}
	void stream() {
	
		int next;
		std::cout << "size=" << size;
		auto ck = dq::c;
#pragma omp parallel for private(next)
		for (int i = 0; i < size; i++) {
			if (c[i].type == 1) {
				for (int k = 1; k <= half; k++) {
					next = i + ck[k] * mask;
					if (next < 0) { continue; }
					if (next >= size) { continue; }

				

					assert(next >= 0);
					assert(next < size);
					//std::cout << c[i].opposite(k) << "swap" << c[next][k] << std::endl;
					std::swap(c[i].opposite(k), c[next][k]);


				}
			}
		}
		
	
	}

	void collide_stream() {
		Vectorxd<d, float> u;
		float rho;
		Vectorxd<q, float> feq;
		auto ck = dq::c;

		//#pragma omp parallel for
		for (int i = 0; i < size; i++) {

			if (c[i].type == 1) {

				phys::collide(c[i]);

				//c[i].swap();
			}
			if (c[i].type == 3) {

				phys::collide(c[i]);

				c[i].swap();
			}
			if (c[i].type == 4) {

				phys::collide(c[i]);

				c[i].swap();
			}
			int next ;
			for (int k = 1; k <= half; k++) {
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
	inline void upwind_swap(int next, int i, int k) {
			float ftmp = c[i][k];
			c[i][k] = c[i].opposite(k);
			c[i].opposite(k) = c[next][k];
			c[next][k] = ftmp;
		}
	void cuInit() {
		//if (hpc->mpi_rank == 0){
			std::cout << "call initCu="<<(*hpc).mpi_rank<<"\n";
			
			sogol::Cell<dq::d, dq::q>* cell=initCu<dq>(c,size);
			delete[] c;
			c = cell;
			for (int i = 0; i < size; i++) {
				std::cout << c[i];
			}
			for (auto b : bcs) {
				b->cell=c;
			}
		

		//}

	}
	void cuRun() {
		//if (hpc->mpi_rank == 0) {
			
			runCu<dq,BGK<DQ<2,9>>>(c, size,mask);
			

		//}

	}


	~Lattice() {
		delete[] c;
		
	}
	

};

template<typename dq, typename phys>
inline std::ostream& operator << (std::ostream& out, Lattice<dq, phys> grid) {
	for (int i = 0; i < grid.size; i++) {
		out << grid.c[i]<<" ";
	};
	return out;
}




}
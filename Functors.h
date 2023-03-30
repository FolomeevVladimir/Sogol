#pragma once
#include <iostream>
#include "Vectorxd.h"
#include <cstdlib>
#include <vector>
#include "Cell.h"
#include <ostream>
#include "Physics.h"
#include <string>

namespace sogol {
	template<unsigned d, unsigned q>
	struct func {
		virtual void operator()(sogol::Cell<d, q> &c, sogol::Vectorxd<d, unsigned> ind) = 0;
	};
	template<unsigned d, unsigned q>
	struct F: func<d,q> {
		virtual void operator()(sogol::Cell<d, q> &c, sogol::Vectorxd<d, unsigned> ind) override {
			std::cout << ind << " "<<c<<std::endl;
		}
	};
	template<unsigned d, unsigned q>
	struct Settype : func<d, q> {
		Vectorxd<d, unsigned> l;
		Settype(Vectorxd<d, unsigned> L) {
			l = L;
			std::cout << l << std::endl;
		}
		virtual void operator()(sogol::Cell<d, q> &c, sogol::Vectorxd<d, unsigned> ind) override {
			
			
			if (ind[1] == l[1] - 1) { c.type = 2; }
			if (ind[1] == 0) { c.type = 2; }
			
			if (ind[0] == 0 && ind[1]!=0 && ind[1]!=(l[1]-1)) { c.type = 3; }
			if (ind[0] == l[0]-1 && ind[1] != 0 && ind[1] != (l[1] - 1)) { c.type = 4; }
		
		}
		
	};
	template<unsigned d, unsigned q>
	struct Save : func<d, q> {
		sogol::BGK<sogol::D2Q9> bgk;
		std::ofstream file;
		Vectorxd<d, unsigned> l;
		
		void init(Vectorxd<d, unsigned> L, int rank) { 
			l = L;
			std::string s = std::to_string(rank)+"u.csv";
			
			
			file.open(s); }

		virtual void operator()(sogol::Cell<d, q>& c, sogol::Vectorxd<d, unsigned> ind) override {
			//if (ind[0]!=0&&ind[1]!=0&&ind[0]!=l[0]-1&& ind[1] != l[1] - 1) {
			if (c.type==1){

				file << ind[0]<< " " <<ind[1]<<" "<<(bgk.u(c))[0]<<" "<< (bgk.u(c))[1]<<" "<<bgk.rho(c)<<std::endl;
			}
		}
	};
}
#pragma once

#include "Vectorxd.h"
#include "Cell.h"
#include "Stencil.h"
//#include "Lattice.h"
#include <vector>
#include "HPC.h"




namespace sogol {

	template<typename Stencil>
	struct BC {
	
		static constexpr unsigned d = Stencil::d;
		static constexpr unsigned q = Stencil::q;
		static constexpr double cs2 = Stencil::cs2;
		Vectorxd<q, double> w = Stencil::w;
		Vectorxd<q, Vectorxd<d, int>> c = Stencil::c;
		int type;
		std::vector<int> indexes;
		Cell<d, q>* cell;
		Vectorxd<d, unsigned> mask;

		BC() {};
		BC(Vectorxd<d, unsigned> mask, Cell<d, q>* c, int type) { this->mask = mask; this->cell = c; this->type=type; }

		inline virtual void process( int i) = 0;
		inline virtual void process() = 0;
	};

	template<typename Stencil>
	struct Periodic :public BC<typename Stencil> {
		static constexpr unsigned d = Stencil::d;
		static constexpr unsigned q = Stencil::q;
		static constexpr unsigned half = (q - 1) / 2;
		static constexpr double cs2 = Stencil::cs2;
		Vectorxd<q, double> w = Stencil::w;
		Vectorxd<q, Vectorxd<d, int>> c = Stencil::c;
		
		
		Vectorxd<d, int> period;

		Periodic(Vectorxd<d, unsigned> mask, Cell<d, q>* c, int type, Vectorxd<d, int> per) : BC<Stencil>::BC(mask,c,type), period(per) { }
		
		inline virtual void process(int i) {
			int next;
		
			auto ck = Stencil::c;
			for (unsigned k = 1; k <= half; k++) {
				
				next = i + ck[k] *(this->mask);

			
					if ((ck[k]*period) < 0) { next = next + period*(this->mask); }
					assert(next >= 0);
		
					std::swap(this->cell[i].opposite(k), this->cell[next][k]);
			
			}
		
		}
		inline virtual void process() {
			for (int i : this->indexes) {
				this->process(i);
			}
		}

	};
	template<typename Stencil>
	struct MPIexchange :public BC<typename Stencil> {
		static constexpr unsigned d = Stencil::d;
		static constexpr unsigned q = Stencil::q;
		static constexpr unsigned half = (q - 1) / 2;
		static constexpr double cs2 = Stencil::cs2;
		Vectorxd<q, double> w = Stencil::w;
		Vectorxd<q, Vectorxd<d, int>> c = Stencil::c;
		sogol::HPC* hpc;

		Vectorxd<d, int> period;

		MPIexchange(Vectorxd<d, unsigned> mask, Cell<d, q>* c, int type, Vectorxd<d, int> per, sogol::HPC* h) : BC<Stencil>::BC(mask, c, type), period(per), hpc(h) { }

		inline virtual void process(int i) {
			int next;
		
			auto ck = Stencil::c;
			for (unsigned k = 1; k <= half; k++) {

				next = i + ck[k] * (this->mask);

				if ((ck[k] * period) < 0) { next = next + period * (this->mask); }
				
				assert(next >= 0);
			
				hpc->swap(&this->cell[i].opposite(k), &this->cell[next][k]);

			}

		}
		inline virtual void process() {
			for (int i : this->indexes) {
				this->process(i);
			}
		}

	};

}
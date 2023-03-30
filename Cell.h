#pragma once
#include <cmath>
#include <iostream>
#include <assert.h> 
#include "Vectorxd.h"
namespace sogol{

template<unsigned d,unsigned q>
struct Cell{

	unsigned type;
	sogol::Vectorxd<q, double> f;
	
	inline Cell() {
		type = 1;
		f=sogol::Vectorxd<q,double>();
	}
	
	inline Cell (std::initializer_list<double> list, unsigned t = 1) {
		
		assert(q == list.size()); 
		type = t;
		f=sogol::Vectorxd<q, double>(list);
	}
	inline double& operator [] (unsigned i){
		assert(i < q);
		return f[i];
	}
	inline double& opposite(unsigned i) {
		
		assert(i >= 0);
		assert(i < q);
		if (i == 0) { return f[i]; }
		if (i <= (q - 1) / 2) { return f[i + (q - 1) / 2]; }
		if (i > (q - 1) / 2) { return f[i - (q - 1) / 2]; }
	}
	inline void swap() {
		for (unsigned i = 1; i <= (q - 1) / 2; i++) {
			std::swap(f[i], f[i + (q - 1) / 2]);
		}
	}
	inline void swap(Cell &c) {
		for (unsigned i = 0; i <q; i++) {
				std::swap(f[i ], c[i]);
		}
	
	}
	

};
template<unsigned d,unsigned q>
inline std::ostream& operator << (std::ostream& out, Cell<d,q> c) {
	out << "type=" << c.type << " " << c.f;
	return out;
}


}//end sogol
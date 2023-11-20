#pragma once
#include <cmath>
#include <iostream>
#include <assert.h> 
#include "Vectorxd.h"
namespace sogol{

template<int d,int q>
struct Cell{

	int type;
	sogol::Vectorxd<q, float> f;
	
	inline Cell() {
		type = 1;
		f=sogol::Vectorxd<q,float>();
	}
	
	inline Cell (std::initializer_list<float> list, int t = 1) {
		
		assert(q == list.size()); 
		type = t;
		f=sogol::Vectorxd<q, float>(list);
	}
	inline constexpr float& operator [] (int i){
		assert(i < q);
		return f[i];
	}
	inline constexpr float& opposite(int i) {
		
		assert(i >= 0);
		assert(i < q);
		if (i == 0) { return f[i]; }
		if (i <= (q - 1) / 2) { return f[i + (q - 1) / 2]; }
		if (i > (q - 1) / 2) { return f[i - (q - 1) / 2]; }
	}
	constexpr inline void swap() {
		for (int i = 1; i <= (q - 1) / 2; i++) {
			float t = f[i + (q - 1) / 2];
			f[i + (q - 1) / 2] = f[i];
			f[i] = t;
			//std::swap(f[i], f[i + (q - 1) / 2]);
		}
	}
	constexpr inline void swap(Cell &c) {
		for (int i = 0; i <q; i++) {
			float t = c[i];
			c[i] = f[i];
			f[i] = t;
			//	std::swap(f[i ], c[i]);
		}
	
	}
	

};
template<int d,int q>
inline std::ostream& operator << (std::ostream& out, Cell<d,q> c) {
	out << "type=" << c.type << " " << c.f;
	return out;
}


}//end sogol
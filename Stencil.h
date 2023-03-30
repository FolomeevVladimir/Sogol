#pragma once
#include "Vectorxd.h"


namespace sogol {

	template<unsigned D, unsigned Q>
	struct Stencil {
		static constexpr unsigned d = D;
		static constexpr unsigned q = Q;
		static double cs2;
		static const Vectorxd<q, Vectorxd<d, int>> c;
		static const Vectorxd<q, double> w;
	};

	struct D2Q9 :Stencil<2, 9> {
		 static double constexpr cs2 = 1. / 3.;


	
	};
	 Vectorxd<D2Q9::q, Vectorxd<D2Q9::d, int>> const D2Q9::c = {
	  { 0, 0},
	  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
	  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
	 };
	 Vectorxd<D2Q9::q, double> const D2Q9::w = { 4. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9 };
}


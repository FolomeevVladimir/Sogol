#pragma once
#include "Vectorxd.h"
//#include <concepts>





namespace sogol {


	template<int D, int Q>
	struct DQ {
		static constexpr int d = D;
		static constexpr int q = Q;
		static constexpr float cs2 = 1. / 3.;
		static constexpr Vectorxd<q, Vectorxd<d, int>> c{
	  { 0, 0},
	  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
	  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
		};
		static constexpr Vectorxd<q, float> w = { 4. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9, 1. / 36, 1. / 9 };
	};

}

#pragma once
#include "Vectorxd.h"
#include "Cell.h"
#include "Stencil.h"




namespace sogol {


template<typename DQ>
struct BGK {
	static constexpr int d = DQ::d;
	static constexpr int q = DQ::q;
	static constexpr float cs2 = DQ::cs2;
	//static constexpr Vectorxd<q, float> w = DQ::w;
	//static constexpr Vectorxd<q, Vectorxd<d, int>> c = DQ::c;
	static constexpr float rho0 = 1.;
	static constexpr float	alpha = 0.01;
	static constexpr float  omega = 1.0 / (3. * alpha + 0.5);
	

	// Исправить ошибку

	inline static constexpr float rho(Cell<d, q> const& cell) {
		
		return cell.f.rho();
		
	}
	inline  static constexpr Vectorxd<d, float> u(Cell<d, q> const& cell) {
		auto c = DQ::c;
		auto m = cell.f *c;
		return m*(1/cell.f.rho());
	}
	inline static constexpr  Vectorxd<q, float> feq(Vectorxd<d, float> const & u, float const & rho) {
		float cs2 = DQ::cs2;
		float cu, u2=u*u, cs4=cs2*cs2;
		auto c = DQ::c;
		auto w = DQ::w;
		Vectorxd<q, float> feq;
		for (int k = 0; k < q; k++) {
			cu = c[k] * u;
			feq[k] = rho * w[k] * (1 + (cu / cs2) + (cu*cu / (2 * cs4))  - (u2 / (2 * cs2)));
			
				}
		return feq;
	}
	inline  Vectorxd<q, float> fout(Cell<d, q> & cell, Vectorxd<q, float> feq) {
		cell.f= omega * feq + (1 - omega) * cell.f;
		return cell.f;
	}
	inline static constexpr void collide(Cell<d, q> & cell) {
		auto c = DQ::c;
		auto w = DQ::w;
		float cs2 = DQ::cs2;
		float rho = cell.f.rho(), cu{}, u2{}, cs4 = cs2 * cs2;
		Vectorxd<d, float> rhou{}, u{};
		Vectorxd<q, float> feq;
		
		for (int i = 0; i<d; i++) {
			
			for (int k=0; k<q; k++) {
				rhou[i] = rhou[i]+cell.f[k] * c[k][i];
			}
		}
		u = (1 / rho) * rhou;
		u2 = u * u;
		for (int k = 0; k < q; k++) {
			cu = c[k] * u;
			feq[k] = rho *w[k] * (1 + (cu / cs2) + (cu * cu / (2 * cs4)) - (u2 / (2 * cs2)));

		}
		for (int k = 0; k < q; k++) {
			cell.f[k] = omega * feq[k] + (1 - omega) * cell.f[k];
			
		}

	
	}
	
	 inline static constexpr void zero(Cell<d, q>& cell) {
		auto w = DQ::w;
		auto c = DQ::c;
		cell.f[0] =c[0][0];
	}
};

 //Vectorxd<BGK<Stencil<2,9>>::q, float> const BGK<Stencil<2,9>>::w = Stencil<2,9>::w;
 //const float BGK<Stencil<2, 9>>::cs2 = Stencil<2,9>::cs2;


}

#pragma once
#include "Vectorxd.h"
#include "Cell.h"
#include "Stencil.h"




namespace sogol {

template<typename Stencil>
struct Dynamics {
	
	static constexpr unsigned d = Stencil::d;
	static constexpr unsigned q = Stencil::q;
	static constexpr double cs2 = Stencil::cs2;
	Vectorxd<q,double> w = Stencil::w;
	Vectorxd<q, Vectorxd<d,int>> c = Stencil::c;

	inline virtual double rho(Cell<d,q> const& cell) = 0;
	inline virtual Vectorxd<d,double> u(Cell<d,q> const& cell) = 0;
	inline virtual Vectorxd<q,double> feq(Vectorxd<d,double> const & u, double const & rho) = 0;
	inline virtual Vectorxd<q, double> fout( Cell<d,q> & cell, Vectorxd<q, double> feq) = 0;
	inline virtual void collide(Cell<d, q>& cell) = 0;
};
template<typename Stencil>
struct BGK : public Dynamics<typename Stencil>{
	static constexpr unsigned d = Stencil::d;
	static constexpr unsigned q = Stencil::q;
	static constexpr double cs2 = Stencil::cs2;
	Vectorxd<q, double> w = Stencil::w;
	Vectorxd<q, Vectorxd<d, int>> c = Stencil::c;
	double rho0 = 1.;
	double	alpha = 0.01;
	double  omega = 1.0 / (3. * alpha + 0.5);

	// Исправить ошибку

	inline  double rho(Cell<d, q> const& cell) {
		return cell.f.rho();
		
	}
	inline  Vectorxd<d, double> u(Cell<d, q> const& cell) {
		auto m = cell.f * c;
		return m*(1/cell.f.rho());
	}
	inline  Vectorxd<q, double> feq(Vectorxd<d, double> const & u, double const & rho) {
		double cu, u2=u*u, cs4=cs2*cs2;
		Vectorxd<q, double> feq;
		for (unsigned k = 0; k < q; k++) {
			cu = c[k] * u;
			feq[k] = rho * w[k] * (1 + (cu / cs2) + (cu*cu / (2 * cs4))  - (u2 / (2 * cs2)));
			
				}
		return feq;
	}
	inline Vectorxd<q, double> fout(Cell<d, q> & cell, Vectorxd<q, double> feq) {
		cell.f= omega * feq + (1 - omega) * cell.f;
		return cell.f;
	}
	inline void collide(Cell<d, q> & cell) {

		double rho = cell.f.rho(), cu,u2,cs4=cs2*cs2;
		Vectorxd<d, double> rhou{}, u{};
		Vectorxd<q, double> feq;
		
		for (unsigned i = 0; i<d; i++) {
			
			for (unsigned k=0; k<q; k++) {
				rhou[i] = rhou[i]+cell.f[k] * c[k][i];
			}
		}
		u = (1 / rho) * rhou;
		u2 = u * u;
		for (unsigned k = 0; k < q; k++) {
			cu = c[k] * u;
			feq[k] = rho * w[k] * (1 + (cu / cs2) + (cu * cu / (2 * cs4)) - (u2 / (2 * cs2)));

		}
		for (unsigned k = 0; k < q; k++) {
			cell.f[k] = omega * feq[k] + (1 - omega) * cell.f[k];
		}

	
	}

};





}

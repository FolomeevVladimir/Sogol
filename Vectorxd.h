
#pragma once
#include <cmath>
#include <iostream>
#include <assert.h> 


namespace sogol {

	template<unsigned d = 2, class T = double>
	struct Vectorxd {

	
		T data[d];
		inline Vectorxd(T value = T{}) {
			assert(d >0);
			for (unsigned i = 0; i < d; i++) { data[i] = value; }

		}
		inline Vectorxd(std::initializer_list<T> list) {
			assert(d == list.size());
			assert(d >0);
			unsigned i = 0;
			for (auto elem : list) {
				data[i] = elem;
				i++;
			}

		}
		inline T& operator [] (const unsigned i) {
			assert(i < d);
			return data[i];
		}
		inline double norm() const {
			double sum = 0;
			for (unsigned i = 0; i < d; i++) {
				sum += data[i] * data[i];
			}
			return sqrt(sum);
		}
		inline T conv() const {
			T res = 1;
			for (unsigned i = 0; i < d; i++) {
				res *= data[i];
			}
			return res;
		}
		inline double norm2() const {
			double sum = 0;
			for (unsigned i = 0; i < d; i++) {
				sum += data[i] * data[i];
			}
			return sum;
		}
		inline T rho() const {
			T sum = {};
			for (unsigned i = 0; i < d; i++) {
				sum += data[i];
			}
			return sum;
		}
		

	};
	template<unsigned d, typename T>
	inline std::ostream& operator << (std::ostream& out, Vectorxd<d, T> vec) {
		out << "{ ";
		for(unsigned i = 0; i < d; i++) {
			out << vec.data[i] << " ";
		}
		out << "}"<<std::endl;
		return out;
	}
	template<unsigned d, typename T,typename W>
	inline Vectorxd<d,decltype(T{}+W{})> operator + (const Vectorxd<d, T> &a,const Vectorxd<d, W> &b) {

		Vectorxd < d, decltype(T{}+W{}) > c;
		for (unsigned i = 0; i < d; i++) {
			c[i] = a.data[i] + b.data[i];
			
		}
		return c;
	
	}
	template<unsigned d, typename T, typename W>
	inline decltype(T{}*W{}) operator * (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		decltype(T{}*W{})  c = {};
		for (unsigned i = 0; i < d; i++) {
			c=c+a.data[i]*b.data[i];

		}
		return c;

	}
	template<unsigned d, typename T, typename W>
	inline Vectorxd < d, decltype(T{}*W{}) > operator * (const T & a, const Vectorxd<d, W>& b) {

		Vectorxd < d, decltype(T{}*W{}) > c;
		for (unsigned i = 0; i < d; i++) {
			c[i] = a*b.data[i];

		}
		return c;

	}
	template<unsigned d, typename T, typename W>
	inline Vectorxd < d, decltype(T{}*W{}) > operator * ( const Vectorxd<d, W>& b, const T& a) {

		Vectorxd < d, decltype(T{}*W{}) > c;
		for (unsigned i = 0; i < d; i++) {
			c[i] = a * b.data[i];

		}
		return c;

	}
	template<unsigned d, typename T, typename W>
	inline Vectorxd<d,bool> operator > (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		Vectorxd<d, bool>  c = {};
		for (unsigned i = 0; i < d; i++) {
			c = a.data[i] > b.data[i];

		}
		return c;

	}
	template<unsigned d, typename T, typename W>
	inline Vectorxd<d, bool> operator < (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		Vectorxd<d, bool>  c = {};
		for (unsigned i = 0; i < d; i++) {
			c = a.data[i] < b.data[i];

		}
		return c;

	}











} //end sogol




#pragma once
#include <cmath>
#include <iostream>
#include <assert.h> 


namespace sogol {

	template<int d = 2, class T = float>
	struct Vectorxd {

	
		T data[d]{};
		inline constexpr Vectorxd(T value = T{}) {
			assert(d >0);
			for (int i = 0; i < d; i++) { data[i] = value; }

		}
		inline constexpr Vectorxd(T* ref ) {
			assert(d > 0);
			for (int i = 0; i < d; i++) { data[i] = ref[i]; }

		}
		inline constexpr Vectorxd(std::initializer_list<T> list) {
			assert(d == list.size());
			assert(d >0);
			int i = 0;
			for (auto elem : list) {
				data[i] = elem;
				i++;
			}

		}
		inline constexpr T& operator [] (const int i) {
			assert(i < d);
			return data[i];
		}
		inline constexpr float norm() const {
			float sum = 0;
			for (int i = 0; i < d; i++) {
				sum += data[i] * data[i];
			}
			return sqrt(sum);
		}
		inline constexpr T conv() const {
			T res = 1;
			for (int i = 0; i < d; i++) {
				res *= data[i];
			}
			return res;
		}
		inline constexpr float norm2() const {
			float sum = 0;
			for (int i = 0; i < d; i++) {
				sum += data[i] * data[i];
			}
			return sum;
		}
		inline constexpr T rho() const {
			T sum = {};
			for (int i = 0; i < d; i++) {
				sum += data[i];
			}
			return sum;
		}
		

	};
	template<int d, typename T>
	inline std::ostream& operator << (std::ostream& out, Vectorxd<d, T> vec) {
		out << "{ ";
		for(int i = 0; i < d; i++) {
			out << vec.data[i] << " ";
		}
		out << "}"<<std::endl;
		return out;
	}
	template<int d, typename T,typename W>
	inline constexpr Vectorxd<d,decltype(T{}+W{})> operator + (const Vectorxd<d, T> &a,const Vectorxd<d, W> &b) {

		Vectorxd < d, decltype(T{}+W{}) > c;
		for (int i = 0; i < d; i++) {
			c[i] = a.data[i] + b.data[i];
			
		}
		return c;
	
	}
	template<int d, typename T, typename W>
	inline constexpr decltype(T{}*W{}) operator * (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		decltype(T{}*W{})  c = {};
		for (int i = 0; i < d; i++) {
			c=c+a.data[i]*b.data[i];

		}
		return c;

	}
	template<int d, typename T, typename W>
	inline constexpr Vectorxd < d, decltype(T{}*W{}) > operator * (const T & a, const Vectorxd<d, W>& b) {

		Vectorxd < d, decltype(T{}*W{}) > c;
		for (int i = 0; i < d; i++) {
			c[i] = a*b.data[i];

		}
		return c;

	}
	template<int d, typename T, typename W>
	inline constexpr Vectorxd < d, decltype(T{}*W{}) > operator * ( const Vectorxd<d, W>& b, const T& a) {

		Vectorxd < d, decltype(T{}*W{}) > c;
		for (int i = 0; i < d; i++) {
			c[i] = a * b.data[i];

		}
		return c;

	}
	template<int d, typename T, typename W>
	inline constexpr Vectorxd<d,bool> operator > (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		Vectorxd<d, bool>  c = {};
		for (int i = 0; i < d; i++) {
			c = a.data[i] > b.data[i];

		}
		return c;

	}
	template<int d, typename T, typename W>
	inline constexpr Vectorxd<d, bool> operator < (const Vectorxd<d, T>& a, const Vectorxd<d, W>& b) {

		Vectorxd<d, bool>  c = {};
		for (int i = 0; i < d; i++) {
			c = a.data[i] < b.data[i];

		}
		return c;

	}











} //end sogol



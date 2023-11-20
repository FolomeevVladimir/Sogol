#pragma once
#include "cuda_runtime.h"
#include "Engine.cuh"
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>

//#include "Stencil.h"
template<typename dq,typename p>
__global__ void kernel(sogol::Cell<dq::d, dq::q>* c,int size) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if ((c[i].type != 2)&&(i<size)) {

		p::collide(c[i]);

		c[i].swap();
		
	}
	__syncthreads();
}
template<typename dq, typename p>
__global__ void kernel2(sogol::Cell<dq::d, dq::q>* c, int size, sogol::Vectorxd<dq::d, int> mask) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if ((c[i].type != 2) && (i < size)) {

		p::collide(c[i]);
		
		c[i].swap();
		

	}
	__syncthreads();
	int next;
	int half = (dq::q - 1) / 2;
	auto ck = dq::c;
		if ((c[i].type == 1)&&(i<size)) {
			for (int k = 1; k <= half; k++) {
			    next = i + ck[k] * mask;
				if (next < 0) { continue; }
				if (next >= size) { continue; }
				__syncthreads();



				assert(next >= 0);
				assert(next < size);
				//swap
				float t = c[next][k];
				c[next][k] =c[i].opposite(k);
				c[i].opposite(k) = t;


			}
		}
	
}
template<typename dq>
sogol::Cell<dq::d, dq::q>*  initCu(sogol::Cell<dq::d, dq::q>* cell, int size) {
	std::cout << "init engine";
	int device = -1;
	int bytes = size * sizeof(sogol::Cell<dq::d, dq::q>);
	cudaGetDevice(&device);
	sogol::Cell<dq::d, dq::q>* c;
	cudaMallocManaged(& c, bytes);
	std::memcpy(c, cell, bytes);
	cell = c;
	cudaMemPrefetchAsync(&cell, size * sizeof(sogol::Cell<dq::d, dq::q>), device);
	cudaDeviceSynchronize();
	std::cout << "engine ready";
	return cell;

}
template sogol::Cell<2,9>* initCu<sogol::DQ<2,9>>(sogol::Cell<2,9>* c, int size);
 
template<typename dq,typename p>
void runCu(sogol::Cell<dq::d, dq::q>* cell, int size, sogol::Vectorxd<dq::d, int> mask) {
	std::cout << "run engine";
	int device = -1;
	int bytes = size * sizeof(sogol::Cell<dq::d, dq::q>);
	cudaGetDevice(&device);
	kernel<dq,p> << <1 + (size / 16), 16 >> > (cell, size);
	cudaDeviceSynchronize();
}
template void runCu<sogol::DQ<2,9>,sogol::BGK<sogol::DQ<2,9>>>(sogol::Cell<2, 9>* c, int size, sogol::Vectorxd<2, int> mask);



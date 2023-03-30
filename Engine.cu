#pragma once
#include "cuda_runtime.h"
#include "Engine.cuh"
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
template<unsigned d, unsigned q>
__global__ void kernel(sogol::Cell<d, q>* c,int size) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < size) { c[i].f.data[0] += 10; }
	__syncthreads();
}
template<unsigned d, unsigned q>
sogol::Cell<d, q>*  initCu(sogol::Cell<d, q>* cell, int size) {
	std::cout << "init engine";
	int device = -1;
	int bytes = size * sizeof(sogol::Cell<d, q>);
	cudaGetDevice(&device);
	sogol::Cell<d, q>* c;
	cudaMallocManaged(& c, bytes);
	std::memcpy(c, cell, bytes);
	cell = c;
	cudaMemPrefetchAsync(&cell, size * sizeof(sogol::Cell<d, q>), device);
	cudaDeviceSynchronize();
	std::cout << "engine ready";
	return cell;

}
template sogol::Cell<2, 9>* initCu<2,9>(sogol::Cell<2, 9>* c, int size);
 
template<unsigned d, unsigned q>
void runCu(sogol::Cell<d, q>* cell, int size) {
	std::cout << "run engine";
	int device = -1;
	int bytes = size * sizeof(sogol::Cell<d, q>);
	cudaGetDevice(&device);
	kernel<d, q> << <1 + (size / 16), 16 >> > (cell, size);
	cudaDeviceSynchronize();
}
template void runCu<2, 9>(sogol::Cell<2, 9>* c, int size);



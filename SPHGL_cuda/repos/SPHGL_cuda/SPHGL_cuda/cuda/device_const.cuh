#ifndef DEVICE_CONST_CUH
#define DEVICE_CONST_CUH

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace DeviceConst {
	__constant__ double PI = 3.14159265358979323846;

	__constant__ int p1 = 73856093;
	__constant__ int p2 = 19349663;
	__constant__ int p3 = 83492791;

	__constant__ float dt = 0.01f;

	__constant__ float g = -9.82f;

	__constant__ float m = 0.02f;
	__constant__ float c = 343.0f;
	__constant__ float p = 101325.0f;
	__constant__ float k = 3.0f;
	__constant__ float mu = 3.5f; //mikro
	__constant__ float rho0 = 998.29f;
	__constant__ float surfTension = 0.0728f;
	__constant__ float threshold = 7.065f;
	__constant__ float cr = 0.2f;
}

#endif
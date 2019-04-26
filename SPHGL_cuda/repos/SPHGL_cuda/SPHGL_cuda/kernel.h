#ifndef KERNEL_H
#define KERNEL_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "thrust\device_vector.h"
#include "device_functions.h"
#include "device_atomic_functions.h"

#include <vector>

#include "cutil_math.h"
#include "sph\Const.hpp"

namespace cuda {
	void initSimulation();
	void simulationStep(float Dt);
	void retrievePositionData(float3* positionBufferCPU);
	void testSpatialQuery(float3* positionBufferCPU, int pn);
}
#endif
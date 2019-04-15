#include "kernel.h"

namespace DeviceConst {
	__constant__ float PI = 3.14159265358979323846;

	__constant__ int p1 = 73856093;
	__constant__ int p2 = 19349663;
	__constant__ int p3 = 83492791;

	__constant__ float g = -9.82f;

	__constant__ float m = 0.02f;
	__constant__ float c = 343.0f;
	__constant__ float p = 101325.0f;
	__constant__ float k = 3.0f;
	__constant__ float mu = 3.5f; //mikro
	__constant__ float rho0 = 998.29f;
	__constant__ float surfTension = 0.0728f;
	__constant__ float threshold = 7.065f;
}

//DEVICE FUNCTIONS

namespace device {
	__device__
		int spatialHash3D(float3 pos, float h, unsigned int nH) {
		long long xor2 = static_cast<long long>(static_cast<long long>(pos.x / h) * DeviceConst::p1)
			^ static_cast<long long>(static_cast<long long>(pos.y / h) * DeviceConst::p2)
			^ static_cast<long long>(static_cast<long long>(pos.z / h) * DeviceConst::p3);

		return (nH + (xor2%nH)) % nH;
	}

	__device__
		float W_def(float3 r, float h) {
		float coef = (315) / (64 * DeviceConst::PI * powf(h, 9));

		float rlen = length(r);

		float weight = h < rlen ? 0.0f : (h*h - (r.x*r.x + r.y*r.y + r.z*r.z));

		return coef * (weight * weight * weight);
	}

	__device__
		float3 W_def_grad(float3 r, float h) {
		float coef = (-945) / (32 * DeviceConst::PI * powf(h, 9));
		float weight = (h*h - (r.x*r.x + r.y*r.y + r.z*r.z));

		return coef * weight * weight * r;
	}

	__device__
		float W_def_lapl(float3 r, float h) {
		float coef = (-945) / (32 * DeviceConst::PI * powf(h, 9));
		float rSQSUM = (r.x*r.x + r.y*r.y + r.z*r.z);

		return coef * (h*h - rSQSUM) * (3 * h*h - 7 * rSQSUM);
	}

	__device__
		float W_press(float3 r, float h) {
		float coef = (15) / (DeviceConst::PI * powf(h, 6));

		float rlen = length(r);

		float weight = rlen > h ? 0.0f : (h - rlen);

		return coef * weight * weight * weight;
	}

	__device__
		float3 W_press_grad(float3 r, float h) {
		float coef = (-45) / (DeviceConst::PI * powf(h, 6));

		float rlen = length(r);

		float weight = (h - rlen);

		return coef * weight * weight * (r / rlen);
	}

	__device__
		float W_visc_lapl(float3 r, float h) {
		float coef = (45) / (DeviceConst::PI * powf(h, 6));

		float rlen = length(r);

		float weight = (h - rlen);

		return coef * weight;
	}

	__device__
		float Gamma(float a, float h) {
		float dp = h / 2;
		return a < (dp) ? (1 - (a / dp)) : 0.0f;
	}

	__device__
		float Xsi(float a, float h) {
		float q = a / h;
		float beta = (0.02 * DeviceConst::c * DeviceConst::c) / a;

		if (0 < q && q < (float)(2 / 3))
			return (float)(2 / 3) * beta;
		else if ((float)(2 / 3) < q && q < 1.0f)
			return beta * (2 * q - 1.5f*q*q);
		else if (1 < q && q < 2)
			return 0.5f * beta * (2 - q) * (2 - q);
		else
			return 0.0f;
	}

	__device__
		float B(float x, float y, float h) {
		return Gamma(y, h) * Xsi(x, h);
	}

	namespace vector {
		__device__
		int push_back(int** arr, int n, int elem) {
			int* new_arr = new int[n + 1];
			//int* new_arr = (int*)malloc((n + 1) * sizeof(int));

			for (int i = 0; i < n; ++i)
				new_arr[i] = (*arr)[i];
			new_arr[n] = elem;

			int* old_arr = *arr;
			*arr = new_arr;

			if (n > 0)
				delete[] old_arr;
			//free(old_arr);
			
			return n + 1;
		}
	}

	__device__
	int spatialQuery(int neighbors[], int idx, float h, float nH, float3* positionBuffer, int* hashElemSizeBuffer, int** hashBuffer) {
		int neighbors_num = 0;

		float3 v_h = make_float3(h, h, h);
		float3 rq = positionBuffer[idx];

		float3 min = make_float3(rq.x - v_h.x, rq.y - v_h.y, rq.z - v_h.z);

		float3 max = make_float3(rq.x + v_h.x, rq.y + v_h.y, rq.z + v_h.z);

		float iter = h / 2;

		//int* hashVals;
		int hashVals_num = 0;
		
		int hashVals[100];
		//int real_size = 0;
		for (float x = min.x; x <= max.x; x += iter) {
			for (float y = min.y; y <= max.y; y += iter) {
				float3 xyz = make_float3(x, y, 0.0f);

				int xyzHash = device::spatialHash3D(xyz, h, nH);

				int contains = 0;
				for (int i = 0; i < hashVals_num; ++i) {
					if (hashVals[i] == xyzHash)
						contains = 1;
				}

				if (contains == 0) {
					//hashVals_num = device::vector::push_back(&hashVals, hashVals_num, xyzHash);
					hashVals[hashVals_num++] = xyzHash;
				}
			}
		}
		
		for (int i = 0; i < hashVals_num; ++i) {
			int hashVal = hashVals[i];

			int* L = hashBuffer[hashVal];

			for (int l = 0; l < hashElemSizeBuffer[hashVal]; ++l) {
				int lidx = L[l];

				if (length(rq - positionBuffer[lidx]) <= h) {
					//neighbors_num = device::vector::push_back(neighbors, neighbors_num, lidx);
					neighbors[neighbors_num++] = lidx;
				}
			}
		}

		return neighbors_num;
	}

	__device__
	float calcMassDensity(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer) {
		float rho = 0.0f;
		float3 r = positionBuffer[idx];

		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			float3 rj = positionBuffer[nidx];	
			rho += DeviceConst::m * W_def(r - rj, h);
		}
		return rho;
	}

	__device__
	float calcPressure(int idx, float* densityBuffer) {
		float rho = densityBuffer[idx];
		return DeviceConst::k * (rho - DeviceConst::rho0);
	}

	__device__
	float3 calcPressureForce(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer, float* pressureBuffer, float* densityBuffer) {
		float rho = densityBuffer[idx];
		float p = pressureBuffer[idx];
		float3 r = positionBuffer[idx];

		float3 Fpressure = make_float3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			if (idx != nidx) {
				float3 rj = positionBuffer[nidx];
				float rhoj = densityBuffer[nidx];
				float pj = pressureBuffer[nidx];

				//Fpressure += ((p / (rho*rho)) + (pj / (rhoj*rhoj))) * Const::particleM * W_press_grad(r - rj, Const::h);
				Fpressure += ((p / (rho*rho)) + (pj / (rhoj*rhoj))) * DeviceConst::m * W_press_grad(r - rj, h);
			}
		}

		return (-rho) * Fpressure;
	}

	__device__
	float3 calcViscosityForce(int idx, float h, int* neighbors, int neighbor_num,
							float3* positionBuffer, float* densityBuffer, float3* velocityBuffer0, float3* velocityBuffer1) {
		float3 vel = (velocityBuffer0[idx] + velocityBuffer1[idx]) / 2.0f;
		float3 r = positionBuffer[idx];

		float3 Fviscosity = make_float3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			if (idx != nidx) {
				float3 rj = positionBuffer[nidx];
				float3 velj = (velocityBuffer0[nidx] + velocityBuffer1[nidx]) / 2.0f;
				float rhoj = densityBuffer[nidx];

				//Fviscosity += (Const::particleM / rhoj) * W_visc_lapl(r - rj, Const::h) * (velj - vel);
				Fviscosity += (DeviceConst::m / rhoj) * W_visc_lapl(r - rj, h) * (velj - vel);
			}
		}

		return DeviceConst::mu * Fviscosity;
	}

	__device__
	float3 calcInternalForces(int idx, float h, int* neighbors, int neighbor_num,
		float3* positionBuffer, float* pressureBuffer, float* densityBuffer, float3* velocityBuffer0, float3* velocityBuffer1)
	{
		float3 Fint = make_float3(0.0f, 0.0f, 0.0f);
		Fint += calcPressureForce(idx, h, neighbors, neighbor_num, positionBuffer, pressureBuffer, densityBuffer);
		Fint += calcViscosityForce(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer, velocityBuffer0, velocityBuffer1);

		return Fint;
	}

	__device__
	float3 calcGravityForce(int idx, float* densityBuffer) {
		return densityBuffer[idx] * make_float3(0.0f, DeviceConst::g, 0.0f);
	}

	__device__
	float3 calcInwardSurfaceNormal(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer, float* densityBuffer) {
		float3 r = positionBuffer[idx];

		float3 normal = make_float3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			float rhoj = densityBuffer[nidx];
			float3 rj = positionBuffer[nidx];

			normal += (DeviceConst::m / rhoj) * W_def_grad(r - rj, h);
		}

		return normal;
	}

	__device__
	float calcSurfaceColorFieldLapl(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer, float* densityBuffer) {
		float3 r = positionBuffer[idx];

		float colorFieldLapl = 0.0f;
		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			float rhoj = densityBuffer[nidx];
			float3 rj = positionBuffer[nidx];

			colorFieldLapl += (DeviceConst::m / rhoj) * W_def_lapl(r - rj, h);
		}

		return colorFieldLapl;
	}

	__device__
	float3 calcTensionForce(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer, float* densityBuffer) {
		float3 n = calcInwardSurfaceNormal(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer);

		float3 Fsurface = make_float3(0.0f, 0.0f, 0.0f);
		if (length(n) >= DeviceConst::threshold) {
			Fsurface = (-DeviceConst::surfTension) * calcSurfaceColorFieldLapl(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer) * (n / length(n));
		}

		return Fsurface;
	}

	__device__
	float3 calcExternalForces(int idx, float h, int* neighbors, int neighbor_num, float3* positionBuffer, float* densityBuffer) {
		float3 Fgravity = calcGravityForce(idx, densityBuffer);
		float3 Fsurface = calcTensionForce(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer);

		return Fgravity + Fsurface;
	}

	
}

//END OF DEVICE FUNCTIONS
//
////KERNEL FUNCTIONS

__global__
void initParticles(const int gridResolution,
				float h,
				float3* positionBuffer,
				float* pressureBuffer,
				float* densityBuffer,
				float3* velocityBuffer0,
				float3* velocityBuffer1,
				float3* forceBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {
		
		int idx = id.x + id.y * gridResolution;

		float x = (h / 3)* (id.x - ((gridResolution - 1) / 2));
		float y = (h / 3) * (id.y - ((gridResolution - 1) / 2));

		positionBuffer[idx] = make_float3(x, y, 0.0f);
		pressureBuffer[idx] = DeviceConst::p;
		densityBuffer[idx] = DeviceConst::rho0;
		velocityBuffer0[idx] = make_float3(0.0f, 0.0f, 0.0f);
		velocityBuffer1[idx] = make_float3(0.0f, 0.0f, 0.0f);
		forceBuffer[idx] = make_float3(0.0f, 0.0f, 0.0f);
	}
}

__global__
void testSpatial(const int gridResolution,
				float h,
				float nH,
				int i,
				float3* positionBuffer,
				int* hashElemSizeBuffer,
				int** hashBuffer,
				float3* spatialTestBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		spatialTestBuffer[idx] = make_float3(1.5f, 1.5f, 0.0f);

		if (idx == i) {
			//thrust::device_vector<int> neighbors = device::spatialQuery(idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);
			int neighbors[100];
			int neighbor_num = 0;

			neighbor_num = device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

			for (int i = 0; i < neighbor_num; ++i) {
				//printf("%d\n", neighbors[i]);
				int nidx = neighbors[i];
				spatialTestBuffer[nidx] = positionBuffer[nidx];
			}
		}
	}
}

__global__
void calcDensityAndPressure(int gridResolution,
					float h,
					float nH,
					float3* positionBuffer,
					float* pressureBuffer,
					float* densityBuffer,
					int* hashElemSizeBuffer,
					int** hashBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		int neighbors[100];
		int neighbor_num = 0;
		neighbor_num = device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

		float rho = device::calcMassDensity(idx, h, neighbors, neighbor_num, positionBuffer);
		float p = device::calcPressure(idx, densityBuffer);

		densityBuffer[idx] = rho;
		pressureBuffer[idx] = p;

		//printf("ok\n");
	}
}

__global__
void calcForces(int gridResolution,
				float h,
				float nH,
				float3* positionBuffer,
				float* pressureBuffer,
				float* densityBuffer,
				float3* velocityBuffer0,
				float3* velocityBuffer1,
				float3* forceBuffer,
				int* hashElemSizeBuffer,
				int** hashBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		int neighbors[100];
		int neighbor_num = 0;
		neighbor_num = device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

		float3 internalF = device::calcInternalForces(idx, h, neighbors, neighbor_num, positionBuffer, pressureBuffer, densityBuffer, velocityBuffer0, velocityBuffer1);
		float3 externalF = device::calcExternalForces(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer);

		forceBuffer[idx] = internalF + externalF;

		//printf("ok\n");
	}
}

//END OF KERNEL FUNCTIONS
//
////BUFFERS

float3* d_positionBuffer;
float3* d_spatialTestBuffer;

float* d_pressureBuffer;

float* d_densityBuffer;

int inputVelocityBuffer = 0;
float3 * d_velocityBuffer[2];

int inputForceBuffer = 0;
float3 * d_forceBuffer[2];

int* d_hashElemSizeBuffer;
int** d_hashBuffer;

//END OF BUFFERS
//
////GENERIC

int gridResolution = Const::partNumX;

dim3 threadsPerBlock(11, 11);
dim3 numBlocks(gridResolution / threadsPerBlock.x, gridResolution / threadsPerBlock.y);

//END OF GENERIC
//
////
void errorCheck();


int spatialHash3D(float3 pos, float h, unsigned int nH) {
	long long xor2 = static_cast<long long>(static_cast<long long>(pos.x / h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>(pos.y / h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>(pos.z / h) * Const::p3);

	return (nH + (xor2%nH)) % nH;
}

void allocateMemory() {
	
	//Position
	cudaMalloc((void**)&d_positionBuffer, sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_spatialTestBuffer, sizeof(float3)*gridResolution*gridResolution);
	
	//Pressure
	cudaMalloc((void**)&d_pressureBuffer, sizeof(float)*gridResolution*gridResolution);

	//Density
	cudaMalloc((void**)&d_densityBuffer, sizeof(float)*gridResolution*gridResolution);

	//Velocity
	cudaMalloc((void**)&d_velocityBuffer[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_velocityBuffer[1], sizeof(float3)*gridResolution*gridResolution);

	//Force
	cudaMalloc((void**)&d_forceBuffer[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_forceBuffer[1], sizeof(float3)*gridResolution*gridResolution);

	//Hash
	cudaMalloc((void**)&d_hashElemSizeBuffer, sizeof(int)*Const::nH);
	cudaMalloc((void**)&d_hashBuffer, sizeof(int*) * Const::nH);
}

void initHashBuffer() {
	float3* positions = new float3[gridResolution*gridResolution];

	cudaMemcpy(positions, d_positionBuffer, sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);

	std::vector< std::vector<int> > hash_table = std::vector< std::vector<int> >( Const::nH );

	//printf("%d", hash_table.size());

	for (int i = 0; i < gridResolution*gridResolution; ++i) {

		float3 pos = positions[i];
		int hashIndex = spatialHash3D(pos, Const::h, Const::nH);

		hash_table[hashIndex].push_back(i);
	}

	int** hashBuffer;
	int* hashElemSizeBuffer = new int[Const::nH];
	hashBuffer = new int*[Const::nH];

	for (int i = 0; i < hash_table.size(); ++i) {
		std::vector<int> neighborList = hash_table[i];

		if (neighborList.size() != 0) {
			hashElemSizeBuffer[i] = neighborList.size();

			//hashBuffer[i] = new int[neighborList.size()];
			int* tempNeighbors = new int[neighborList.size()];

			for (int j = 0; j < neighborList.size(); ++j) {
				tempNeighbors[j] = neighborList[j];
			}

			int* d_tempNeighbors;
			cudaMalloc((void**)&d_tempNeighbors, sizeof(int) * neighborList.size());
			cudaMemcpy(d_tempNeighbors, tempNeighbors, sizeof(int) * neighborList.size(), cudaMemcpyHostToDevice);

			hashBuffer[i] = d_tempNeighbors;

			delete[] tempNeighbors;
		}
		else {
			hashElemSizeBuffer[i] = 0;
		}

		//printf("%d\n", hashElemSizeBuffer[i]);
	}

	//printf("%d", hashElemSizeBuffer[1287]);

	cudaMemcpy(d_hashElemSizeBuffer, hashElemSizeBuffer, sizeof(int)*Const::nH, cudaMemcpyHostToDevice);
	cudaMemcpy(d_hashBuffer, hashBuffer, sizeof(int*) * Const::nH, cudaMemcpyHostToDevice);
}

void cuda::initSimulation() {
	allocateMemory();

	initParticles<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		d_positionBuffer,
		d_pressureBuffer,
		d_densityBuffer,
		d_velocityBuffer[0],
		d_velocityBuffer[1],
		d_forceBuffer[0]
		);

	initHashBuffer();
}

void densityAndPressureStep() {
	calcDensityAndPressure<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		d_positionBuffer,
		d_pressureBuffer,
		d_densityBuffer,
		d_hashElemSizeBuffer,
		d_hashBuffer
		);

	errorCheck();
}

void forcesStep() {
	calcForces<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		d_positionBuffer,
		d_pressureBuffer,
		d_densityBuffer,
		d_velocityBuffer[0],
		d_velocityBuffer[1],
		d_forceBuffer[0],
		d_hashElemSizeBuffer,
		d_hashBuffer
		);

	errorCheck();
}


void cuda::simulationStep() {
	densityAndPressureStep();
	forcesStep();
}

void cuda::retrievePositionData(float3* positionBufferCPU) {
	cudaMemcpy(positionBufferCPU, d_positionBuffer, sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);
}

void cuda::testSpatialQuery(float3* positionBufferCPU, int pn) {
	testSpatial<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		pn,
		d_positionBuffer,
		d_hashElemSizeBuffer,
		d_hashBuffer,
		d_spatialTestBuffer
		);

	errorCheck();
	cudaMemcpy(positionBufferCPU, d_spatialTestBuffer, sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);
}

void errorCheck() {
	// Check for any errors launching the kernel
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Launch failed: %s\n", cudaGetErrorString(cudaStatus));
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d!\n", cudaStatus);
	}
}


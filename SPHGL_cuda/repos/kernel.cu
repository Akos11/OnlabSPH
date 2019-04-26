#include "kernel.h"

namespace DeviceConst {
	__constant__ double PI = 3.14159265358979323846;

	__constant__ int p1 = 73856093;
	__constant__ int p2 = 19349663;
	__constant__ int p3 = 83492791;

	__constant__ float dt = 0.002f;

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

//DEVICE FUNCTIONS

namespace device {

	__device__
	int sgn(float val) {
		return (0.0f < val) - (val < 0.0f);
	}

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

				if (contains == 0 && hashVals_num < 100) {
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

				if (length(rq - positionBuffer[lidx]) <= h && neighbors_num < 100) {
					//neighbors_num = device::vector::push_back(neighbors, neighbors_num, lidx);
					neighbors[neighbors_num++] = lidx;
				}
			}
		}

		return neighbors_num;
	}

	__device__
	int spatialQueryBorder(int neighbors[], int idx, float h, float borderNH, float3 pos, float3* boundaryPositionBuffer, int* boundaryHashElemSizeBuffer, int** boundaryHashBuffer) {
		int neighbors_num = 0;

		float3 v_h = make_float3(h, h, h);
		float3 rq = pos;

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

				int xyzHash = device::spatialHash3D(xyz, h, borderNH);

				int contains = 0;
				for (int i = 0; i < hashVals_num; ++i) {
					if (hashVals[i] == xyzHash)
						contains = 1;
				}

				if (contains == 0 && hashVals_num < 100) {
					//hashVals_num = device::vector::push_back(&hashVals, hashVals_num, xyzHash);
					hashVals[hashVals_num++] = xyzHash;
				}
			}
		}

		for (int i = 0; i < hashVals_num; ++i) {
			int hashVal = hashVals[i];

			int* L = boundaryHashBuffer[hashVal];

			for (int l = 0; l < boundaryHashElemSizeBuffer[hashVal]; ++l) {
				int lidx = L[l];
				float3 boundaryPos = boundaryPositionBuffer[lidx];
				
				if (idx == 234)
					;// printf("%f %f\n", boundaryPos.x, boundaryPos.y);

				if (length(rq - boundaryPos) <= h && neighbors_num < 100) {
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
							float3* positionBuffer, float* densityBuffer, float3* velocityBufferPrev, float3* velocityBufferCurr) {
		float3 vel = (velocityBufferPrev[idx] + velocityBufferCurr[idx]) / 2.0f;
		float3 r = positionBuffer[idx];

		float3 Fviscosity = make_float3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < neighbor_num; ++i) {
			int nidx = neighbors[i];
			if (idx != nidx) {
				float3 rj = positionBuffer[nidx];
				float3 velj = (velocityBufferPrev[nidx] + velocityBufferCurr[nidx]) / 2.0f;
				float rhoj = densityBuffer[nidx];

				//Fviscosity += (Const::particleM / rhoj) * W_visc_lapl(r - rj, Const::h) * (velj - vel);
				Fviscosity += (DeviceConst::m / rhoj) * W_visc_lapl(r - rj, h) * (velj - vel);
			}
		}

		return DeviceConst::mu * Fviscosity;
	}

	__device__
	float3 calcInternalForces(int idx, float h, int* neighbors, int neighbor_num,
		float3* positionBuffer, float* pressureBuffer, float* densityBuffer, float3* velocityBufferPrev, float3* velocityBufferCurr)
	{
		float3 Fint = make_float3(0.0f, 0.0f, 0.0f);
		Fint += calcPressureForce(idx, h, neighbors, neighbor_num, positionBuffer, pressureBuffer, densityBuffer);
		Fint += calcViscosityForce(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer, velocityBufferPrev, velocityBufferCurr);

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

	__device__
	float calcBoundaryF(float3 x) {
		float3 c = make_float3(0.0f, 0.0f, 0.0f);
		float r = 0.8f;

		float xMinusc = length(x - c);
		return xMinusc*xMinusc - r*r;
	}

	__device__
	float3 getBoundaryContactPoint(float3 x) {
		float3 c = make_float3(0.0f, 0.0f, 0.0f);
		float r = 0.8f;

		return c + r * ((x - c) / length(x - c));
	}

	__device__
	float getBoundaryDepth(float3 x) {
		float3 c = make_float3(0.0f, 0.0f, 0.0f);
		float r = 0.8f;

		return fabsf(length(c - x) - r);
	}

	__device__
	float3 getBoundarySurfaceNormal(float3 x) {
		float3 c = make_float3(0.0f, 0.0f, 0.0f);

		float fx = calcBoundaryF(x);

		return sgn(fx) * ((c - x) / length(c - x));
	}

	__device__
	float3 getVelAfterCollision(float3 vel, float3 n, float depth) {
		float coef = 1 + DeviceConst::cr * (depth / (DeviceConst::dt * length(vel)));
		float udotn = dot(vel, n);

		return vel - coef * udotn * n;
	}
	
}

//END OF DEVICE FUNCTIONS
//
////KERNEL FUNCTIONS

__global__
void initParticles(const int gridResolution,
				float h,
				float3* positionBufferIn,
				float3* positionBufferOut,
				float* pressureBuffer,
				float* densityBufferIn,
				float* densityBufferOut,
				float3* velocityBufferPrevIn,
				float3* velocityBufferPrevOut,
				float3* velocityBufferCurrIn,
				float3* velocityBufferCurrOut,
				float3* forceBufferIn,
				float3* forceBufferOut)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {
		
		int idx = id.x + id.y * gridResolution;

		float x = (h / 3)* (id.x - ((gridResolution - 1) / 2));
		float y = (h / 3) * (id.y - ((gridResolution - 1) / 2));

		positionBufferIn[idx] = make_float3(x, y, 0.0f);
		positionBufferOut[idx] = make_float3(x, y, 0.0f);
		pressureBuffer[idx] = DeviceConst::p;
		densityBufferIn[idx] = DeviceConst::rho0;
		densityBufferOut[idx] = DeviceConst::rho0;
		velocityBufferPrevIn[idx] = make_float3(0.0f, 0.0f, 0.0f);
		velocityBufferPrevOut[idx] = make_float3(0.0f, 0.0f, 0.0f);
		velocityBufferCurrIn[idx] = make_float3(0.0f, 0.0f, 0.0f);
		velocityBufferCurrOut[idx] = make_float3(0.0f, 0.0f, 0.0f);
		forceBufferIn[idx] = make_float3(0.0f, 0.0f, 0.0f);
		forceBufferOut[idx] = make_float3(0.0f, 0.0f, 0.0f);
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
					float* densityBufferIn,
					float* densityBufferOut,
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
		float p = device::calcPressure(idx, densityBufferIn);

		densityBufferOut[idx] = rho;
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
				float3* velocityBufferPrev,
				float3* velocityBufferCurr,
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

		float3 internalF = device::calcInternalForces(idx, h, neighbors, neighbor_num, positionBuffer, pressureBuffer, densityBuffer, velocityBufferPrev, velocityBufferCurr);
		float3 externalF = device::calcExternalForces(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer);

		forceBuffer[idx] = internalF + externalF;

		//printf("ok\n");
	}
}

__global__
void stepTimeIntegrator(int gridResolution,
						float3* positionBufferIn,
						float3* positionBufferOut,
						float* densityBuffer,
						float3* velocityBufferPrevIn,
						float3* velocityBufferPrevOut,
						float3* velocityBufferCurrIn,
						float3* velocityBufferCurrOut,
						float3* forceBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		float3 prevPos = positionBufferIn[idx];
		float3 prevVel = velocityBufferCurrIn[idx];

		velocityBufferPrevOut[idx] = prevVel;

		float3 F = forceBuffer[idx];
		//printf("%f %f\n", F.x, F.y);
		float rho = densityBuffer[idx];

		float3 currVel = prevVel + DeviceConst::dt * (F / rho);
		float3 currPos = prevPos + DeviceConst::dt * currVel;

		velocityBufferCurrOut[idx] = currVel;
		positionBufferOut[idx] = currPos;

		/*if (idx == 135)
			printf("%f %f\n", prevVel.x, prevVel.y);*/
	}
}

__global__
void handleCollisions(int gridResolution,
					float3* positionBufferIn,
					float3* positionBufferOut,
					float3* velocityBufferPrevIn,
					float3* velocityBufferPrevOut,
					float3* velocityBufferCurrIn,
					float3* velocityBufferCurrOut)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		float3 pos = positionBufferIn[idx];

		float boundaryF = device::calcBoundaryF(pos);

		if (boundaryF > 0.0f) {
			float3 contactPoint = device::getBoundaryContactPoint(pos);
			float3 velBeforeCollision = velocityBufferCurrIn[idx];

			float3 surfaceNormal = device::getBoundarySurfaceNormal(contactPoint);
			float depth = device::getBoundaryDepth(pos);	//Sztem ide ez kell es nem a contactPoint

			float3 velAfterCollision = device::getVelAfterCollision(velBeforeCollision, surfaceNormal, depth);

			velocityBufferPrevOut[idx] = velocityBufferPrevIn[idx];
			velocityBufferCurrOut[idx] = (velBeforeCollision + velAfterCollision) / 2.0f;
			positionBufferOut[idx] = contactPoint -0.001f * (contactPoint / length(contactPoint));
		} else {
			velocityBufferPrevOut[idx] = velocityBufferPrevIn[idx];
			velocityBufferCurrOut[idx] = velocityBufferCurrIn[idx];
			positionBufferOut[idx] = positionBufferIn[idx];
		}
	}
}

__global__
void calcBoundaryForces(int gridResolution,
						float h,
						float boundaryNH,
						float3* positionBuffer,
						float3* boundaryPositionBuffer,
						float3* forceBufferIn,
						float3* forceBufferOut,
						int* boundaryHashSizeElemBuffer,
						int** boundaryHashBuffer)
{
	int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

	if (id.x >= 0 && id.x < gridResolution &&
		id.y >= 0 && id.y < gridResolution) {

		int idx = id.x + id.y * gridResolution;

		float3 pPos = positionBuffer[idx];

		int borderNeighbors[100];
		int bneighbor_num = 0;
		bneighbor_num = device::spatialQueryBorder(borderNeighbors, idx, h, boundaryNH, pPos, boundaryPositionBuffer, boundaryHashSizeElemBuffer, boundaryHashBuffer);

		float3 Fak = make_float3(0.0f, 0.0f, 0.0f);

		for (int i = 0; i < bneighbor_num; ++i) {
			int kidx = borderNeighbors[i];
			float3 kPos = boundaryPositionBuffer[kidx];

			//printf("%f %f\n", kPos.x, kPos.y);

			float3 xak = kPos - pPos;
			float3 surfaceNormal = device::getBoundarySurfaceNormal(pPos);

			float y = dot(xak, surfaceNormal);
			float x = sqrtf(length(xak) * length(xak) - y * y);

			Fak += -(DeviceConst::m / (DeviceConst::m + DeviceConst::m)) * device::B(x, y, h) * surfaceNormal;
		}

		float3 inForce = forceBufferIn[idx];

		forceBufferOut[idx] = inForce + Fak;
	}
}

//END OF KERNEL FUNCTIONS
//
////BUFFERS

int inputPositionBuffer = 0;
float3* d_positionBuffer[2];

float3* d_spatialTestBuffer;

float* d_pressureBuffer;

int inputDensityBuffer = 0;
float* d_densityBuffer[2];

int inputVelocityBuffer = 0;
float3 * d_velocityBufferPrev[2];
float3 * d_velocityBufferCurr[2];

int inputForceBuffer = 0;
float3 * d_forceBuffer[2];

int* d_hashElemSizeBuffer;
int** d_hashBuffer;

//Boundary
float3* d_boundaryPositionBuffer;

int* d_boundaryHashElemSizeBuffer;
int** d_boundaryHashBuffer;

int start = 1;

//END OF BUFFERS
//
////GENERIC

int gridResolution = Const::partNumX;

dim3 threadsPerBlock(3, 3);
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

int spatialHash3DBorder(float3 pos) {
	long long xor2 = static_cast<long long>(static_cast<long long>(pos.x / Const::h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>(pos.y / Const::h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>(pos.z / Const::h) * Const::p3);

	return (Const::borderNH + (xor2%Const::borderNH)) % Const::borderNH;
}

void allocateMemory() {
	
	//Position
	cudaMalloc((void**)&d_positionBuffer[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_positionBuffer[1], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_spatialTestBuffer, sizeof(float3)*gridResolution*gridResolution);
	
	//Pressure
	cudaMalloc((void**)&d_pressureBuffer, sizeof(float)*gridResolution*gridResolution);

	//Density
	cudaMalloc((void**)&d_densityBuffer[0], sizeof(float)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_densityBuffer[1], sizeof(float)*gridResolution*gridResolution);

	//Velocity
	cudaMalloc((void**)&d_velocityBufferPrev[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_velocityBufferPrev[1], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_velocityBufferCurr[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_velocityBufferCurr[1], sizeof(float3)*gridResolution*gridResolution);

	//Force
	cudaMalloc((void**)&d_forceBuffer[0], sizeof(float3)*gridResolution*gridResolution);
	cudaMalloc((void**)&d_forceBuffer[1], sizeof(float3)*gridResolution*gridResolution);

	//Hash
	cudaMalloc((void**)&d_hashElemSizeBuffer, sizeof(int)*Const::nH);
	cudaMalloc((void**)&d_hashBuffer, sizeof(int*) * Const::nH);

	////Boundary
	cudaMalloc((void**)&d_boundaryPositionBuffer, sizeof(float3)*Const::borderParticleNum);
	
	cudaMalloc((void**)&d_boundaryHashElemSizeBuffer, sizeof(int)*Const::borderNH);
	cudaMalloc((void**)&d_boundaryHashBuffer, sizeof(int*) * Const::borderNH);
}

int** hashBufferPrev;
int* hashBufferElemSizePrev;

void freePrev() {
	for (int i = 0; i < Const::nH; ++i) {
		if (hashBufferElemSizePrev[i] > 0) {
			cudaFree(hashBufferPrev[i]);
		}
	}
}

void updateHashBuffer() {
	float3* positions = new float3[gridResolution*gridResolution];

	cudaMemcpy(positions, d_positionBuffer[inputPositionBuffer], sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);

	std::vector< std::vector<int> > hash_table = std::vector< std::vector<int> >(Const::nH);

	//printf("%d", hash_table.size());

	for (int i = 0; i < gridResolution*gridResolution; ++i) {

		float3 pos = positions[i];
		int hashIndex = spatialHash3D(pos, Const::h, Const::nH);

		hash_table[hashIndex].push_back(i);
	}

	delete[] positions;

	int* hashElemSizeBuffer = new int[Const::nH];
	int** hashBuffer = new int*[Const::nH];

	if (!start) {
		freePrev();
	}

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

	if (!start) {
		delete[] hashBufferElemSizePrev;
		delete[] hashBufferPrev;
	}

	hashBufferElemSizePrev = hashElemSizeBuffer;
	hashBufferPrev = hashBuffer;
}

void initBoundaryParticles() {
	float3* positions = new float3[Const::borderParticleNum];

	float radStep = (2 * Const::PI) / Const::borderParticleNum;

	for (int i = 0; i < Const::borderParticleNum; ++i) {
		float x = Const::borderR * cosf(radStep * i);
		float y = Const::borderR * sinf(radStep * i);

		//printf("%f %f\n", x, y);

		positions[i] = make_float3(x, y, 0.0f);
	}

	cudaMemcpy(d_boundaryPositionBuffer, positions, sizeof(float3)*Const::borderParticleNum, cudaMemcpyHostToDevice);

	std::vector< std::vector<int> > hash_table = std::vector< std::vector<int> >(Const::borderNH);

	for (int i = 0; i < Const::borderParticleNum; ++i) {

		float3 pos = positions[i];
		int hashIndex = spatialHash3D(pos, Const::h, Const::borderNH);

		hash_table[hashIndex].push_back(i);
	}

	int** boundaryHashBuffer;
	int* boundaryHashElemSizeBuffer = new int[Const::borderNH];
	boundaryHashBuffer = new int*[Const::borderNH];

	for (int i = 0; i < hash_table.size(); ++i) {
		std::vector<int> neighborList = hash_table[i];

		if (neighborList.size() != 0) {
			boundaryHashElemSizeBuffer[i] = neighborList.size();

			//printf("%d\n", neighborList.size());

			int* tempNeighbors = new int[neighborList.size()];

			for (int j = 0; j < neighborList.size(); ++j) {
				tempNeighbors[j] = neighborList[j];
			}

			int* d_tempNeighbors;
			cudaMalloc((void**)&d_tempNeighbors, sizeof(int) * neighborList.size());
			cudaMemcpy(d_tempNeighbors, tempNeighbors, sizeof(int) * neighborList.size(), cudaMemcpyHostToDevice);

			boundaryHashBuffer[i] = d_tempNeighbors;

			delete[] tempNeighbors;
		}
		else {
			boundaryHashElemSizeBuffer[i] = 0;
		}
	}

	cudaMemcpy(d_boundaryHashElemSizeBuffer, boundaryHashElemSizeBuffer, sizeof(int)*Const::borderNH, cudaMemcpyHostToDevice);
	cudaMemcpy(d_boundaryHashBuffer, boundaryHashBuffer, sizeof(int*) * Const::borderNH, cudaMemcpyHostToDevice);
}

void cuda::initSimulation() {
	allocateMemory();

	initParticles<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		d_positionBuffer[0],
		d_positionBuffer[1],
		d_pressureBuffer,
		d_densityBuffer[0],
		d_densityBuffer[1],
		d_velocityBufferPrev[0],
		d_velocityBufferPrev[1],
		d_velocityBufferCurr[0],
		d_velocityBufferCurr[1],
		d_forceBuffer[0],
		d_forceBuffer[1]
		);

	initBoundaryParticles();

	updateHashBuffer();
	start = 0;

	errorCheck();
}

void densityAndPressureStep() {
		calcDensityAndPressure <<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		d_positionBuffer[inputPositionBuffer],
		d_pressureBuffer,
		d_densityBuffer[inputDensityBuffer],
		d_densityBuffer[(inputDensityBuffer + 1) % 2],
		d_hashElemSizeBuffer,
		d_hashBuffer
		);

		inputDensityBuffer = (inputDensityBuffer + 1) % 2;

	errorCheck();
}

void forcesStep() {
	calcForces<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		d_positionBuffer[inputPositionBuffer],
		d_pressureBuffer,
		d_densityBuffer[inputDensityBuffer],
		d_velocityBufferPrev[inputVelocityBuffer],
		d_velocityBufferCurr[inputVelocityBuffer],
		d_forceBuffer[(inputForceBuffer + 1) % 2],
		d_hashElemSizeBuffer,
		d_hashBuffer
		);

	inputForceBuffer = (inputForceBuffer + 1) % 2;

	errorCheck();
}

void timeStep() {
	stepTimeIntegrator<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		d_positionBuffer[inputPositionBuffer],
		d_positionBuffer[(inputPositionBuffer + 1) % 2],
		d_densityBuffer[inputDensityBuffer],
		d_velocityBufferPrev[inputVelocityBuffer],
		d_velocityBufferPrev[(inputVelocityBuffer + 1) % 2],
		d_velocityBufferCurr[inputVelocityBuffer],
		d_velocityBufferCurr[(inputVelocityBuffer + 1) % 2],
		d_forceBuffer[inputForceBuffer]
		);

	inputPositionBuffer = (inputPositionBuffer + 1) % 2;
	inputVelocityBuffer = (inputVelocityBuffer + 1) % 2;

	errorCheck();
}

void collisionHandlingStep() {
	handleCollisions<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		d_positionBuffer[inputPositionBuffer],
		d_positionBuffer[(inputPositionBuffer + 1) % 2],
		d_velocityBufferPrev[inputVelocityBuffer],
		d_velocityBufferPrev[(inputVelocityBuffer + 1) % 2],
		d_velocityBufferCurr[inputVelocityBuffer],
		d_velocityBufferCurr[(inputVelocityBuffer + 1) % 2]
		);

	inputPositionBuffer = (inputPositionBuffer + 1) % 2;
	inputVelocityBuffer = (inputVelocityBuffer + 1) % 2;

	errorCheck();
}

void boundaryForcesStep() {
	calcBoundaryForces<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::borderNH,
		d_positionBuffer[inputPositionBuffer],
		d_boundaryPositionBuffer,
		d_forceBuffer[inputForceBuffer],
		d_forceBuffer[(inputForceBuffer + 1) % 2],
		d_boundaryHashElemSizeBuffer,
		d_boundaryHashBuffer
		);

	inputForceBuffer = (inputForceBuffer + 1) % 2;

	errorCheck();
}

void cuda::simulationStep() {
	densityAndPressureStep();
	forcesStep();
	boundaryForcesStep();
	timeStep();

	collisionHandlingStep();

	updateHashBuffer();
}

void cuda::retrievePositionData(float3* positionBufferCPU) {
	cudaMemcpy(positionBufferCPU, d_positionBuffer[inputPositionBuffer], sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);
}

void cuda::testSpatialQuery(float3* positionBufferCPU, int pn) {
	testSpatial<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Const::h,
		Const::nH,
		pn,
		d_positionBuffer[inputPositionBuffer],
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



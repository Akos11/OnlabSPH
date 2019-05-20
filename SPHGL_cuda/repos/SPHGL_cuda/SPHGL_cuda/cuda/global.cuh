#ifndef GLOBAL_CUH
#define GLOBAL_CUH

#include "device.cuh"

namespace Global {

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
			int tempidX = id.x % 33;
			int tempidY = id.y % 33;
			float x = (h / 3)* (tempidX - ((gridResolution/33 - 1) / 2))* 1.2f;
			float y = (h / 3) * (tempidY - ((gridResolution/33 - 1) / 2)) * 1.2f;
			//négyzetbe inicializálás
			/*
			float x = (h / 3)* (id.x - ((gridResolution - 1) / 2));
			float y = (h / 3) * (id.y - ((gridResolution - 1) / 2));
			*/

			float z = (h / 3) * ((id.x / 33) + (id.y / 33) * (gridResolution / 33))*3.5f;// +0.2f;


			positionBufferIn[idx] = make_float3(x, y, z);
			positionBufferOut[idx] = make_float3(x, y, z);
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
				//thrust::Device_vector<int> neighbors = Device::spatialQuery(idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);
				int neighbors[100];
				int neighbor_num = 0;

				neighbor_num = Device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

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
			neighbor_num = Device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

			float rho = Device::calcMassDensity(idx, h, neighbors, neighbor_num, positionBuffer);
			float p = Device::calcPressure(idx, rho);
			//float p = 0.0f;
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
			neighbor_num = Device::spatialQuery(neighbors, idx, h, nH, positionBuffer, hashElemSizeBuffer, hashBuffer);

			float3 internalF = Device::calcInternalForces(idx, h, neighbors, neighbor_num, positionBuffer, pressureBuffer, densityBuffer, velocityBufferPrev, velocityBufferCurr);
			float3 externalF = Device::calcExternalForces(idx, h, neighbors, neighbor_num, positionBuffer, densityBuffer);

			forceBuffer[idx] = internalF + externalF;

			/*if (idx == 324)
				printf("%f\n", length(forceBuffer[idx]));*/

				//printf("ok\n");
		}
	}

	__global__
		void stepTimeIntegrator(int gridResolution,
			float Dt,
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

			float3 currVel = prevVel + Dt * (F / rho);
			float3 currPos = prevPos + Dt * currVel;

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

			float boundaryF = Device::calcBoundaryF(pos);

			if (boundaryF > 0.0f) {
				float3 contactPoint = Device::getBoundaryContactPoint(pos);
				float3 velBeforeCollision = velocityBufferCurrIn[idx];

				float3 surfaceNormal = Device::getBoundarySurfaceNormal(contactPoint);
				float depth = Device::getBoundaryDepth(pos);	//Sztem ide ez kell es nem a contactPoint

				float3 velAfterCollision = Device::getVelAfterCollision(velBeforeCollision, surfaceNormal, depth);

				velocityBufferPrevOut[idx] = velocityBufferPrevIn[idx];
				velocityBufferCurrOut[idx] = (velBeforeCollision + velAfterCollision) / 2.0f;
				positionBufferOut[idx] = contactPoint - 0.001f * (contactPoint / length(contactPoint));
			}
			else {
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
			bneighbor_num = Device::spatialQueryBorder(borderNeighbors, idx, h, boundaryNH, pPos, boundaryPositionBuffer, boundaryHashSizeElemBuffer, boundaryHashBuffer);

			float3 Fak = make_float3(0.0f, 0.0f, 0.0f);

			for (int i = 0; i < bneighbor_num; ++i) {
				int kidx = borderNeighbors[i];
				float3 kPos = boundaryPositionBuffer[kidx];

				//printf("%f %f\n", kPos.x, kPos.y);

				float3 xak = kPos - pPos;
				float3 surfaceNormal = Device::getBoundarySurfaceNormal(kPos);

				float y = dot(xak, surfaceNormal);
				float x = sqrtf(length(xak) * length(xak) - y * y);

				Fak += -(DeviceConst::m / (DeviceConst::m + DeviceConst::m)) * Device::B(x, y, h) * surfaceNormal;
			}

			float3 inForce = forceBufferIn[idx];

			forceBufferOut[idx] = inForce + 0.01f * Fak;
		}
	}

	__global__
		void updateHash(int gridResolution,
			float h,
			unsigned int nH,
			float3* positionBuffer,
			int* hashElemSizeBuffer,
			int** hashBuffer) {

		for (int i = 0; i < nH; ++i) {
			hashElemSizeBuffer[i] = 0;
		}

		for (int i = 0; i < gridResolution*gridResolution; ++i) {
			float3 pos = positionBuffer[i];

			int hash = Device::spatialHash3D(pos, h, nH);
			int elemSize = hashElemSizeBuffer[hash];
			if (elemSize < 100) {
				hashBuffer[hash][elemSize] = i;
				hashElemSizeBuffer[hash]++;
			}
		}
	}

	__global__
		void updateHashParallel(int gridResolution,
			float h,
			unsigned int nH,
			float3* positionBuffer,
			int* hashElemSizeBuffer,
			int** hashBuffer)
	{
		int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

		if (id.x >= 0 && id.x < gridResolution &&
			id.y >= 0 && id.y < gridResolution) {

			int idx = id.x + id.y * gridResolution;

			float3 pos = positionBuffer[idx];
			int hash = Device::spatialHash3D(pos, h, nH);


			int nextIdx = atomicAdd(&hashElemSizeBuffer[hash], 1);
			//printf("%d\n", nextIdx);
			if (nextIdx < 100) {
				hashBuffer[hash][nextIdx] = idx;
			}
			else {
				atomicSub(&hashElemSizeBuffer[hash], 1);
			}
		}
	}

	__global__
		void resetHashElemSizeBuffer(int gridResolution, int lastIdx, int nH, int* hashElemSizeBuffer) {
		int2 id = make_int2(blockIdx.x * blockDim.x + threadIdx.x, blockIdx.y * blockDim.y + threadIdx.y);

		if (id.x >= 0 && id.x < gridResolution &&
			id.y >= 0 && id.y < gridResolution) {

			int idx = id.x + id.y * gridResolution;

			if (idx < lastIdx) {
				hashElemSizeBuffer[idx] = 0;
			}
			else {
				int restNum = nH - lastIdx;

				for (int i = 0; i < restNum; ++i) {
					hashElemSizeBuffer[lastIdx + i] = 0;
				}
			}
		}
	}

}
//END OF KERNEL FUNCTIONS

#endif

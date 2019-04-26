#ifndef DEVICE_CUH
#define DEVICE_CUH

#include "..\cutil_math.h"
#include "device_const.cuh"

namespace Device {

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
		float coef = (315.0f) / (64 * DeviceConst::PI * powf(h, 9));

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

				int xyzHash = Device::spatialHash3D(xyz, h, nH);

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

				int xyzHash = Device::spatialHash3D(xyz, h, borderNH);

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
		float calcPressure(int idx, float rho) {
		//float rho = densityBuffer[idx];
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
		return xMinusc * xMinusc - r * r;
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

#endif
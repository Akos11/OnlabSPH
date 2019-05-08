#include "kernel.h"

#include "cuda\global.cuh"

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

	const int nH = Const::nH;
	int** hashLists = new int*[Const::nH];
	for (int i = 0; i < Const::nH; ++i) {
		int* d_hashList;
		cudaMalloc((void**)&d_hashList, sizeof(int) * 100);
		hashLists[i] = d_hashList;
	}

	cudaMemcpy(d_hashBuffer, hashLists, sizeof(int*) * Const::nH, cudaMemcpyHostToDevice);

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

void update() {
	Global::updateHash<<<1, 1>>>(gridResolution,
						Const::h,
						Const::nH,
						d_positionBuffer[inputPositionBuffer],
						d_hashElemSizeBuffer,
						d_hashBuffer);
}

int maxPowerOfTwo = static_cast<int>(std::log2(Const::nH));
int b = static_cast<int>(sqrtf(pow(2, maxPowerOfTwo) / 64.0f));
int lastIdx = b * b * 64 - 1;

void updateParallel() {
	dim3 tpb(8, 8);
	//printf("%lf\n", pow(2, maxPowerOfTwo));
	dim3 nb(b, b);

	
	Global::resetHashElemSizeBuffer<<<nb, tpb>>>(8 * b, lastIdx, Const::nH, d_hashElemSizeBuffer);
	Global::updateHashParallel<<<numBlocks, threadsPerBlock>>>(gridResolution,
														Const::h,
														Const::nH,
														d_positionBuffer[inputPositionBuffer],
														d_hashElemSizeBuffer,
														d_hashBuffer);

	errorCheck();
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

void initBoundaryParticles3D() {
	float3* positions = new float3[Const::borderParticleNum];

	float radStep = (2 * Const::PI) / Const::borderParticleNum;

	int bordNum = 0;
	float hHalf = Const::h / 2;
	int thetaNum = (2 * Const::PI) / (Const::h / 2);
	int fiNum = (Const::PI) / (Const::h / 2);
	for (int t = 0; t < thetaNum; ++t) {
		for (int f = 0; f < fiNum; ++f) {
			float x = Const::borderR * cosf(t * hHalf) * sinf(f * hHalf);
			float y = Const::borderR * sinf(t * hHalf) * sinf(f * hHalf);
			float z = Const::borderR * cosf(f * hHalf);

			positions[bordNum] = make_float3(x, y, z);
			
			bordNum++;
		}
	}

	cudaMemcpy(d_boundaryPositionBuffer, positions, sizeof(float3)*Const::borderParticleNum, cudaMemcpyHostToDevice);

	std::vector< std::vector<int> > hash_table = std::vector< std::vector<int> >(Const::borderNH);

	for (int i = 0; i < Const::borderParticleNum; ++i) {

		float3 pos = positions[i];
		int hashIndex = spatialHash3D(pos, Const::h, Const::borderNH);

		hash_table[hashIndex].push_back(i);
	}

	int** boundaryHashBuffer = new int*[Const::borderNH];
	int* boundaryHashElemSizeBuffer = new int[Const::borderNH];

	for (int i = 0; i < hash_table.size(); ++i) {
		std::vector<int> neighborList = hash_table[i];

		if (neighborList.size() != 0) {
			boundaryHashElemSizeBuffer[i] = neighborList.size();

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

	errorCheck();
}

void cuda::initSimulation() {
	allocateMemory();

	Global::initParticles<<<numBlocks, threadsPerBlock>>>(
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

	initBoundaryParticles3D();

	updateParallel();
	//update();
	//updateHashBuffer();
	start = 0;

	errorCheck();
}

void densityAndPressureStep() {
	Global::calcDensityAndPressure <<<numBlocks, threadsPerBlock>>>(
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
	Global::calcForces<<<numBlocks, threadsPerBlock>>>(
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

void timeStep(float Dt) {
	Global::stepTimeIntegrator<<<numBlocks, threadsPerBlock>>>(
		gridResolution,
		Dt,
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
	Global::handleCollisions<<<numBlocks, threadsPerBlock>>>(
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
	Global::calcBoundaryForces<<<numBlocks, threadsPerBlock>>>(
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

void cuda::simulationStep(float Dt) {
	densityAndPressureStep();
	forcesStep();
	boundaryForcesStep();
	timeStep(Dt);
	
	collisionHandlingStep();

	//updateHashBuffer();
	//update();
	updateParallel();
}

void cuda::retrievePositionData(float3* positionBufferCPU) {
	cudaMemcpy(positionBufferCPU, d_positionBuffer[inputPositionBuffer], sizeof(float3) * gridResolution * gridResolution, cudaMemcpyDeviceToHost);
}

void cuda::testSpatialQuery(float3* positionBufferCPU, int pn) {
	Global::testSpatial<<<numBlocks, threadsPerBlock>>>(
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



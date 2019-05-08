#ifndef CONST_HPP
#define CONST_HPP

#include "Vec.hpp"
#include <cmath>

static bool isPrime(int n) {
	if (n <= 3)
		return n > 1;
	else if (((n % 2) == 0) || ((n % 3) == 0))
		return false;

	int i = 5;
	while (i * i <= n) {
		if (((n%i) == 0) || ((n % (i + 2)) == 0))
			return false;
		i = i + 6;
	}
	return true;
}

static unsigned int genNextPrime(int n) {
	unsigned int num = n;
	bool isprime = isPrime(num);
	while (!isprime) {
		num++;
		isprime = isPrime(num);
	}

	return num;
}

static unsigned int getSphereParticleNum(float R, float h, double PI) {
	int thetaNum = (2 * PI) / (h / 2);
	int fiNum = (PI) / (h / 2);

	return thetaNum * fiNum;
}

namespace Const {
	const bool DDD = true;
	const bool box = false;

	const double PI = 3.14159265358979323846;

	const Vec3 g = Vec3{ 0.0f, -9.82f, 0.0f };
	const float dt = 0.01f;
	const float p = 101325.0f;
	const float c = 343;

	const float rho0 = 998.29f;
	const float particleM = 0.02f;
	const float mu = 3.5f;	//mikro
	const float surfTension = 0.0728f;
	const float threshold = 7.065f;
	const float cr = 0.0f;
	const float k = 3.0f;
	const float x = 20.0f;
	//	const float x = 10.0f;

	//	const float h = 0.0460f;
	//	const float h = 0.07597f;

	const int p1 = 73856093;
	const int p2 = 19349663;
	const int p3 = 83492791;

	//	const unsigned int particleNum = 4900; //At kell irni a particles initet is
	//	const unsigned int nH = 9803;

	//	const unsigned int nH = 883;	//21*21
	//	const unsigned int nH = 4051;	//45*45

	const float borderR = 0.8f;


	const float gridMin = -0.1f;
	const float gridMax = 0.1f;

	const int partNumX = 60;
	const int partNumY = 60;

	const unsigned int particleNum = partNumX * partNumY;

	const unsigned int nH = genNextPrime(2 * particleNum);//2179;	//33*33

	const float h = cbrtf((3 * 0.1*x) / (4 * PI*particleNum));

	//	const unsigned int borderParticleNum = (2 * borderR*PI) / (h / 3);
	const unsigned int borderParticleNum = DDD ? getSphereParticleNum(borderR, h, PI) : (box ? (4 * 1.4) / (h / 3) : (2 * borderR*PI) / (h / 3));
	const unsigned int borderNH = genNextPrime(2 * borderParticleNum);		//263, 397, 541
	//const float gridMin = -1;
	//const float gridMax = 1;
}
#endif

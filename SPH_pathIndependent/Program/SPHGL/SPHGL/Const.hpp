#ifndef CONST_HPP
#define CONST_HPP

#include "Vec.hpp"
#include <cmath>

namespace Const {
	const bool DDD = true;

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

	const unsigned int nH = 2179;	//33*33
//	const unsigned int nH = 883;	//21*21
//	const unsigned int nH = 4051;	//45*45

	const float borderR = 0.8f;


	const float gridMin = -0.1f;
	const float gridMax = 0.1f;

	const int partNumX = 33;
	const int partNumY = 33;

	const unsigned int particleNum = partNumX * partNumY;

	const float h = cbrtf((3 * 0.1*x) / (4 * PI*particleNum));

	const unsigned int borderParticleNum = (2 * borderR*PI) / (h / 3);
	const unsigned int borderNH = 397;		//263, 397, 541
	//const float gridMin = -1;
	//const float gridMax = 1;
}
#endif

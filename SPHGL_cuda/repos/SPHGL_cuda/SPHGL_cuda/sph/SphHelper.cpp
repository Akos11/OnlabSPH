//#include "stdafx.h"
#include "SphHelper.hpp"

float randFloatBtw(float min, float max) {
	float rnd = min + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (max - min)));
	return rnd;
}

int spatialHash3D(const Particle * p) {
	long long xor2 = static_cast<long long>(static_cast<long long>((p->pos).x / Const::h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>((p->pos).y / Const::h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>((p->pos).z / Const::h) * Const::p3);

	return (Const::nH + (xor2%Const::nH)) % Const::nH;
}

int spatialHash3D(const Vec3 r) {
	long long xor2 = static_cast<long long>(static_cast<long long>(r.x / Const::h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>(r.y / Const::h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>(r.z / Const::h) * Const::p3);

	return (Const::nH + (xor2%Const::nH)) % Const::nH;
}

int spatialHash3DBorder(const Particle * p) {
	long long xor2 = static_cast<long long>(static_cast<long long>((p->pos).x / Const::h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>((p->pos).y / Const::h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>((p->pos).z / Const::h) * Const::p3);

	return (Const::borderNH + (xor2%Const::borderNH)) % Const::borderNH;
}

int spatialHash3DBorder(const Vec3 r) {
	long long xor2 = static_cast<long long>(static_cast<long long>(r.x / Const::h) * Const::p1)
		^ static_cast<long long>(static_cast<long long>(r.y / Const::h) * Const::p2)
		^ static_cast<long long>(static_cast<long long>(r.z / Const::h) * Const::p3);

	return (Const::borderNH + (xor2%Const::borderNH)) % Const::borderNH;
}

////Particle

unsigned int Particle::nextID = 0;

Vec3 Particle::getAcceleration() {
	return this->F / this->rho;
}

Vec3 Particle::getVelocityAtT() const {
	return (this->prevVel + this->currVel) / 2.0f;
}

////END OF Particle

/////Particles

void Particles::initRnd() {

	for (int i = 0; i < Const::particleNum; ++i) {
		float x = randFloatBtw(Const::gridMin, Const::gridMax);
		float y = randFloatBtw(Const::gridMin, Const::gridMax);

		Particle * p = new Particle{ Vec3{ x, y, 0.0f } };
		particles[i] = p;
		insertParticle(p);
	}
}

void Particles::init() {
	for (int i = 0; i < Const::partNumX; ++i) {
		for (int j = 0; j < Const::partNumY; ++j) {
			float x = (Const::h / 3)* (i - ((Const::partNumX - 1) / 2));
			float y = (Const::h / 3) * (j - ((Const::partNumY - 1) / 2));

			Particle * p = new Particle{ Vec3{ x + 0.0f, y, 0.0f } };
			particles[i * Const::partNumX + j] = p;
			insertParticle(p);
		}
	}

	std::cout << "nH: " << Const::nH << std::endl;
}

void Particles::initBorder() {
	if (!Const::DDD) {
		float radStep = (2 * Const::PI) / Const::borderParticleNum;

		int num;
		for (int i = 0; i < Const::borderParticleNum; ++i) {
			float x = Const::borderR * cosf(radStep * i);
			float y = Const::borderR * sinf(radStep * i);

			Particle * p = new Particle{ Vec3{ x, y, 0.0f } };
			borderParticles[i] = p;
			insertBorderParticle(p);
			//num = i;
		}

		std::cout << "Border Particle Num: " << Const::borderParticleNum << " borderNh: " << Const::borderNH << std::endl;
	}
	else {

		int bordNum = 0;
		float hHalf = Const::h / 2;
		int thetaNum = (2 * Const::PI) / (Const::h / 2);
		int fiNum = (Const::PI) / (Const::h / 2);
		for (int t = 0; t < thetaNum; ++t) {
			for (int f = 0; f < fiNum; ++f) {
				float x = Const::borderR * cosf(t * hHalf) * sinf(f * hHalf);
				float y = Const::borderR * sinf(t * hHalf) * sinf(f * hHalf);
				float z = Const::borderR * cosf(f * hHalf);

				Particle * p = new Particle{ Vec3{ x, y, z } };
				borderParticles[bordNum] = p;
				insertBorderParticle(p);
				bordNum++;
			}
		}

		std::cout << bordNum;
	}

	/*
	for (int i = 0; i < borderHash_table.size(); ++i) {
	std::cout << borderHash_table[i].size() << std::endl;
	}
	*/
}

void Particles::initBorderBox() {
	float w = 1.4f;
	float h = 1.4f;

	int numPerSide = Const::borderParticleNum / 4;

	float step = Const::h / 3;

	int n = 0;
	for (int i = 0; i < 4; ++i) {
		if (i == 0) {
			for (int j = 0; j < numPerSide; ++j) {
				float x = w / 2;
				float y = (j * step) - (h / 2);

				Particle * p = new Particle{ Vec3{ x + 0.0f, y, 0.0f } };
				borderParticles[n] = p;
				insertBorderParticle(p);
				n++;
			}

		}
		else if (i == 1) {
			for (int j = 0; j < numPerSide; ++j) {
				float x = -w / 2;
				float y = (j * step) - (h / 2);

				Particle * p = new Particle{ Vec3{ x + 0.0f, y, 0.0f } };
				borderParticles[n] = p;
				insertBorderParticle(p);
				n++;
			}

		}
		if (i == 2) {
			for (int j = 0; j < numPerSide; ++j) {
				float y = h / 2;
				float x = (j * step) - (w / 2);

				Particle * p = new Particle{ Vec3{ x + 0.0f, y, 0.0f } };
				borderParticles[n] = p;
				insertBorderParticle(p);
				n++;
			}

		}
		if (i == 3) {
			for (int j = 0; j < numPerSide; ++j) {
				float y = -h / 2;
				float x = -((j * step) - (w / 2));

				Particle * p = new Particle{ Vec3{ x + 0.0f, y, 0.0f } };
				borderParticles[n] = p;
				insertBorderParticle(p);
				n++;
			}

		}
	}
	borderParticles[n - 1] = new Particle{ Vec3{ -0.68f, -0.7f, 0.0f } };
	insertBorderParticle(borderParticles[n - 1]);
	borderParticles[n] = new Particle{ Vec3{ -0.67f, -0.7f, 0.0f } };
	insertBorderParticle(borderParticles[n]);

	std::cout << n << " " << Const::borderParticleNum;
}

void Particles::insertParticle(Particle * p) {
	int hashIndex = spatialHash3D(p);

	hash_table[hashIndex].push_back(p);
}

void Particles::insertBorderParticle(Particle * p) {
	int hashIndex = spatialHash3DBorder(p);

	borderHash_table[hashIndex].push_back(p);
}

std::vector<Particle *> Particles::spatialQuery(Particle * queryP) {
	std::vector<Particle *> neighbors = std::vector<Particle *>{};

	Vec3 v_h = Vec3{ Const::h, Const::h, Const::h };
	Vec3 rq = Vec3{ queryP->pos };

	float xMin = rq.x - v_h.x;
	float yMin = rq.y - v_h.y;
	float zMin = rq.z - v_h.z;

	float xMax = rq.x + v_h.x;
	float yMax = rq.y + v_h.y;
	float zMax = rq.z + v_h.z;

	float iter = Const::h / 2;

	std::set<int> hashVals = std::set<int>{};
	for (float x = xMin; x <= xMax; x += iter) {
		for (float y = yMin; y <= yMax; y += iter) {
			if (!Const::DDD) {
				Vec3 xyz = Vec3{ x, y, 0.0f };

				hashVals.insert(spatialHash3D(xyz));
			}
			else {
				for (float z = zMin; z <= zMax; z += iter) {
					Vec3 xyz = Vec3{ x, y, z };

					hashVals.insert(spatialHash3D(xyz));
				}
			}
		}
	}

	for (int hashVal : hashVals) {
		std::vector<Particle *> L = hash_table[hashVal];

		for (auto p : L) {
			if ((queryP->pos - p->pos).len() <= Const::h)
				neighbors.push_back(p);
		}
	}

	return neighbors;
}

std::vector<Particle *> Particles::spatialQueryBorder(Particle * queryP) {
	std::vector<Particle *> neighbors = std::vector<Particle *>{};

	Vec3 v_h = Vec3{ Const::h, Const::h, Const::h };
	Vec3 rq = Vec3{ queryP->pos };

	float xMin = rq.x - v_h.x;
	float yMin = rq.y - v_h.y;
	float zMin = rq.z - v_h.z;

	float xMax = rq.x + v_h.x;
	float yMax = rq.y + v_h.y;
	float zMax = rq.z + v_h.z;

	float iter = Const::h / 2;

	std::set<int> hashVals = std::set<int>{};
	for (float x = xMin; x <= xMax; x += iter) {
		for (float y = yMin; y <= yMax; y += iter) {
			if (!Const::DDD) {
				Vec3 xyz = Vec3{ x, y, 0.0f };

				hashVals.insert(spatialHash3DBorder(xyz));
			}
			else {
				for (float z = zMin; z <= zMax; z += iter) {
					Vec3 xyz = Vec3{ x, y, z };

					hashVals.insert(spatialHash3DBorder(xyz));
				}
			}
		}
	}

	for (int hashVal : hashVals) {
		std::vector<Particle *> L = borderHash_table[hashVal];

		for (auto p : L) {
			if ((queryP->pos - p->pos).len() <= 1 * Const::h) {
				neighbors.push_back(p);
			}
		}
	}

	return neighbors;
}

void Particles::updateSpatialHashing() {
	for (int i = 0; i < Const::nH; ++i) {
		hash_table[i].clear();
	}

	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles[i];

		insertParticle(p);
	}

}


////END OF Particles


////Box

Vec3 Box::xToLocal(const Vec3& x) {
	return x - c;
}

Vec3 Box::getContactPoint(const Vec3& x) {
	Vec3 maxExtXoc = callFuncForComponents(std::fmaxf, -ext, x);

	return callFuncForComponents(std::fminf, ext, maxExtXoc);
}

float Box::F(const Vec3& x) {
	Vec3 X = x;
	Vec3 xMinusExt = X.callFuncForComponents(std::fabsf) - ext;

	return xMinusExt.getMaxComponent();
}

float Box::getDepth(const Vec3& x) {
	return (getContactPoint(x) - x).len();
}

Vec3 Box::getSurfaceNormal(const Vec3& x) {
	/*Vec3 notUnitSurfNorm = (getContactPoint(x) - x).callFuncForComponents(sgn<float>);

	float len = notUnitSurfNorm.len();

	if (len == 0)
	len = 0.00001f;

	return notUnitSurfNorm / len; */

	Vec3 norm = Vec3{};
	if (fabs(x.x) > ext.x - 0.01)
		norm.x = -sgn<float>(x.x);
	if (fabs(x.y) > ext.y - 0.01)
		norm.y = -sgn<float>(x.y);

	if (norm.len() > 1)
		norm = norm / norm.len();

	return norm;
}

Vec3 Box::velAfterCollision(const Vec3& vel, const Vec3& n, float depth) {
	float coef = 1 + Const::cr * (depth / (Const::dt*vel.len()));
	float udotn = vel.dot(n);

	return vel - 1 * udotn * n;
}

Vec3 Box::getNormal(const Vec3& x) {
	Vec3 norm = Vec3{};
	if (fabs(x.x) > ext.x - 0.01)
		norm.x = sgn<float>(x.x);
	if (fabs(x.y) > ext.y - 0.01)
		norm.y = sgn<float>(x.y);

	if (norm.len() > 1)
		norm = norm / norm.len();
	return norm;
}
////END OF Box

////Sphere

float Sphere::F(const Vec3& x) {
	float xMinusc = (x - c).len();
	return xMinusc*xMinusc - r*r;
}

Vec3 Sphere::getContactPoint(const Vec3& x) {
	return c + r * ((x - c) / (x - c).len());
}

float Sphere::getDepth(const Vec3& x) {
	return fabsf((c - x).len() - r);
}

Vec3 Sphere::getSurfaceNormal(const Vec3& x) {
	float fx = F(x);

	//if (fx == 0)
	//	fx = 1;

	return sgn(fx) * ((c - x) / (c - x).len());
}

Vec3 Sphere::velAfterCollision(const Vec3& vel, const Vec3& n, float depth) {
	float coef = 1 + Const::cr * (depth / (Const::dt*vel.len()));
	float udotn = vel.dot(n);

	return vel - 1.0f * udotn * n;
}

Vec3 Sphere::getNormal(const Vec3& x) {
	return (x - c) / (x - c).len();
}

////END OF Sphere

////LeapFrogIntegrator

void LeapFrogIntegrator::stepIntegrator(std::vector<Particle *>& particles) {

	//	std::cout << "leap" << std::endl;

	this->t += dt;

	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles[i];

		Vec3 prevPos = p->pos;	//At t - dt
		Vec3 prevVel = p->currVel;	//At t - 1/2dt

		p->prevVel = prevVel;
		Vec3 currVel = prevVel + dt * p->getAcceleration();	//At t + 1/2dt
		Vec3 currPos = prevPos + dt * currVel;				//At t

		p->currVel = currVel;
		p->pos = currPos;
	}
}

////END OF LeapFrogIntegrator

float W_def(const Vec3& r, float h) {
	float coef = (315) / (64 * Const::PI * powf(h, 9));

	//float weight = h < r.len() ? 0.0f : (h*h - (r.len()*r.len()));
	float weight = h < r.len() ? 0.0f : (h*h - (r.x*r.x + r.y*r.y + r.z*r.z));

	return coef * (weight * weight * weight);
}

Vec3 W_def_grad(const Vec3& r, float h) {
	float coef = (-945) / (32 * Const::PI * powf(h, 9));
	float weight = (h*h - (r.x*r.x + r.y*r.y + r.z*r.z));

	return coef * weight * weight * r;
}

float W_def_lapl(const Vec3& r, float h) {
	float coef = (-945) / (32 * Const::PI * powf(h, 9));
	float rSQSUM = (r.x*r.x + r.y*r.y + r.z*r.z);

	return coef * (h*h - rSQSUM) * (3 * h*h - 7 * rSQSUM);
}

float W_press(const Vec3& r, float h) {
	float coef = (15) / (Const::PI * powf(h, 6));

	float weight = r.len() > h ? 0.0f : (h - r.len());

	return coef * weight * weight * weight;
}

Vec3 W_press_grad(const Vec3& r, float h) {
	float coef = (-45) / (Const::PI * powf(h, 6));

	float weight = (h - r.len());

	return coef * weight * weight * (r / r.len());
}

float W_visc_lapl(const Vec3& r, float h) {
	float coef = (45) / (Const::PI * powf(h, 6));

	float weight = (h - r.len());

	return coef * weight;
}

float Gamma(float a) {
	float dp = Const::h / 2;
	return a < (dp) ? (1 - (a / dp)) : 0.0f;
}

float Xsi(float a) {
	float q = a / Const::h;
	float beta = (0.02 * Const::c * Const::c) / a;

	if (0 < q && q < (float)(2 / 3))
		return (float)(2 / 3) * beta;
	else if ((float)(2 / 3) < q && q < 1.0f)
		return beta * (2 * q - 1.5f*q*q);
	else if (1 < q && q < 2)
		0.5f * beta * (2 - q) * (2 - q);
	else
		return 0.0f;
}

float B(float x, float y) {
	return Gamma(y) * Xsi(x);
}



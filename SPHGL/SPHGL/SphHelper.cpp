#include "stdafx.h"
#include "SphHelper.hpp"


float randFloatBtw(float min, float max) {
	float rnd = min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
	return rnd;
}

int spatialHash3D(const Particle * p) {
	long long xor = static_cast<long long>(static_cast<long long>((p->pos).x / Const::h) * Const::p1)
						^ static_cast<long long>(static_cast<long long>((p->pos).y / Const::h) * Const::p2)
						^ static_cast<long long>(static_cast<long long>((p->pos).z / Const::h) * Const::p3);
	
	return (Const::nH + (xor%Const::nH)) % Const::nH;
}

int spatialHash3D(const Vec3 r) {
	long long xor = static_cast<long long>(static_cast<long long>(r.x / Const::h) * Const::p1)
						^ static_cast<long long>(static_cast<long long>(r.y / Const::h) * Const::p2)
						^ static_cast<long long>(static_cast<long long>(r.z / Const::h) * Const::p3);
	
	return (Const::nH + (xor%Const::nH)) % Const::nH;
}

////Particle

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

		Particle * p = new Particle{ Vec3{x, y, 0.0f} };
		particles[i] = p;
		insertParticle(p);
	}
}

void Particles::init() {
	for (int i = 0; i < Const::partNumX; ++i) {
		for (int j = 0; j < Const::partNumY; ++j) {
			float x = (Const::h / 3)* (i-((Const::partNumX-1)/2));
			float y = (Const::h / 3) * (j-((Const::partNumY-1)/2));

			Particle * p = new Particle{ Vec3{x, y, 0.0f} };
			particles[i * Const::partNumX + j] = p;
			insertParticle(p);
		}
	}
}

void Particles::insertParticle(Particle * p) {
	int hashIndex = spatialHash3D(p);

	hash_table[hashIndex].push_back(p);
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
	for (float x = xMin; x <= xMax; x+=iter) {
		for (float y = yMin; y <= yMax; y+=iter) {
				Vec3 xyz = Vec3{ x, y, 0.0f };

				hashVals.insert(spatialHash3D(xyz));
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

Vec3 Box::xToLocal(Vec3 x) {
	return x - c;
}

Vec3 Box::getContactPoint(Vec3 x) {
	Vec3 maxExtXoc = callFuncForComponents(std::fmaxf, -ext, xToLocal(x));

	return callFuncForComponents(std::fminf, ext, maxExtXoc) + c;
}

float Box::F(Vec3 x) {
	Vec3 xLocMinusExt = xToLocal(x).callFuncForComponents(std::fabsf) - ext;

	return xLocMinusExt.getMaxComponent();
}

float Box::getDepth(Vec3 x) {
	return (getContactPoint(x) - x).len();
}

Vec3 Box::getSurfaceNormal(Vec3 x) {
	Vec3 notUnitSurfNorm = ((getContactPoint(x) - c) - xToLocal(x)).callFuncForComponents(sgn<float>);

	float len = notUnitSurfNorm.len();

	if (len == 0)
		len = 0.00001f;

	return notUnitSurfNorm / len;
}

Vec3 Box::velAfterCollision(const Vec3& vel, const Vec3& n, float depth) {
	float coef = 1 + Const::cr * (depth / (Const::dt*vel.len()));
	float udotn = vel.dot(n);
	
	return vel - coef * udotn * n;
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

	return sgn(-fx) * ((c - x) / (c - x).len());
}

Vec3 Sphere::velAfterCollision(const Vec3& vel, const Vec3& n, float depth) {
	float coef = 1 + Const::cr * (depth / (Const::dt*vel.len()));
	float udotn = vel.dot(n);
	
	return vel - coef * udotn * n;
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

	return coef * (h*h - rSQSUM) * (3*h*h - 7*rSQSUM);
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


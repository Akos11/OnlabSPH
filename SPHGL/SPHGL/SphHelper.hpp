#ifndef SPHHELPER_HPP
#define SPHHELPER_HPP

#include "Vec.hpp"
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include "Const.hpp"

float randFloatBtw(float min, float max);

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Particle {

	Vec3 pos;
	Vec3 currVel;	//At t + 1/2dt
	Vec3 prevVel;	//At t - 1/2dt
	Vec3 F;
	float p;
	float rho;
	Vec3 Fint = Vec3{};
	Vec3 Fext = Vec3{};

	Particle(Vec3 pos, Vec3 vel, Vec3 F) : pos{pos}, currVel{vel}, F{F}, p{ Const::p }, rho{ Const::rho0 } {
		currVel = currVel - 0.5f * Const::dt * (F / rho);
	}
	Particle(Vec3 pos) : pos{ pos }, currVel{ Vec3{0.0f, 0.0f, 0.0f} }, F{ Vec3{0.0f, 0.0f, 0.0f} }, p{ Const::p }, rho{ Const::rho0 } {
		currVel = currVel - 0.5f * Const::dt * (F / rho);
	}

	Vec3 getAcceleration();
	Vec3 getVelocityAtT() const;
};

int spatialHash3D(const Particle * p);
int spatialHash3D(const Vec3 r);

struct Particles {

	std::vector<Particle *> particles = std::vector<Particle *>{ Const::particleNum };
	std::vector<std::vector<Particle *> > hash_table = std::vector<std::vector<Particle *> >{ Const::nH };

	Particles() {}

	void initRnd();
	void init(); 

	void insertParticle(Particle * p);

	std::vector<Particle *> spatialQuery(Particle * queryP);

	void updateSpatialHashing();
};

struct Box {
	Vec3 c;
	float width; //x
	float height; //y
	float depth; //z
	Vec3 ext;

	Box(Vec3 c = Vec3{}, float width = 0.9f, float height = 0.9f, float depth = 0)
		: c{ c }, width{ width }, height{ height }, depth{ depth } {
		ext = Vec3{ width, height, depth };
	}

	Vec3 xToLocal(Vec3 x);
	Vec3 getContactPoint(Vec3 x);
	float F(Vec3 x);
	float getDepth(Vec3 x);
	Vec3 getSurfaceNormal(Vec3 x);
	Vec3 velAfterCollision(const Vec3& vel, const Vec3& n, float depth);
};

struct Sphere {
	Vec3 c;
	float r;

	Sphere(Vec3 c = Vec3{}, float r = 0.8f) : c {c}, r{r} {}

	float F(const Vec3& x);
	Vec3 getContactPoint(const Vec3& x);
	float getDepth(const Vec3& x);
	Vec3 getSurfaceNormal(const Vec3& x);
	Vec3 velAfterCollision(const Vec3& vel, const Vec3& n, float depth);
};

struct LeapFrogIntegrator {
	float dt;
	float t;

	LeapFrogIntegrator(float dt = Const::dt, float t = 0.0f) : dt{ dt }, t{ t } {}

	void stepIntegrator(std::vector<Particle *>& particles);
};

float W_def(const Vec3& r, float h);
Vec3 W_def_grad(const Vec3& r, float h);
float W_def_lapl(const Vec3& r, float h);

Vec3 W_press_grad(const Vec3& r, float h);

float W_visc_lapl(const Vec3& r, float h);


#endif // !SPHHELPER_HPP



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
	Vec3 color;

	static unsigned int nextID;
	unsigned int ID	;

	Particle(Vec3 pos, Vec3 vel, Vec3 F) : pos{ pos }, currVel{ vel }, F{ F }, p{ Const::p }, rho{ Const::rho0 } {
		currVel = currVel - 0.5f * Const::dt * (F / rho);
		color = Vec3(randFloatBtw(0.0f, 1.0f), randFloatBtw(0.0f, 1.0f), randFloatBtw(0.0f, 1.0f));
		ID = ++nextID;
	}
	Particle(Vec3 pos) : pos{ pos }, currVel{ Vec3{ 0.0f, 0.0f, 0.0f } }, F{ Vec3{ 0.0f, 0.0f, 0.0f } }, p{ Const::p }, rho{ Const::rho0 } {
		currVel = currVel - 0.5f * Const::dt * (F / rho);
		color = Vec3(randFloatBtw(0.0f, 1.0f), randFloatBtw(0.0f, 1.0f), randFloatBtw(0.0f, 1.0f));
		ID = ++nextID;
	}

	Vec3 getAcceleration();
	Vec3 getVelocityAtT() const;
};

int spatialHash3D(const Particle * p);
int spatialHash3D(const Vec3 r);

int spatialHash3DBorder(const Particle * p);
int spatialHash3DBorder(const Vec3 r);

struct Particles {

	std::vector<Particle *> particles = std::vector<Particle *>{ Const::particleNum };
	std::vector<std::vector<Particle *> > hash_table = std::vector<std::vector<Particle *> >{ Const::nH };

	std::vector<Particle *> borderParticles = std::vector<Particle *>{ Const::DDD ? 13530 : Const::borderParticleNum };
	std::vector<std::vector<Particle *> > borderHash_table = std::vector<std::vector<Particle *> >{ Const::DDD ? 27061 : Const::borderNH };

	Particles() {}

	void initRnd();
	void init();
	void initBorder();
	void initBorderBox();

	void insertParticle(Particle * p);
	void insertBorderParticle(Particle * p);

	std::vector<Particle *> spatialQuery(Particle * queryP);
	std::vector<Particle *> spatialQueryBorder(Particle * queryP);

	void updateSpatialHashing();
};

struct Boundary {
	virtual float F(const Vec3& x) = 0;
	virtual Vec3 getContactPoint(const Vec3& x) = 0;
	virtual float getDepth(const Vec3& x) = 0;
	virtual Vec3 getSurfaceNormal(const Vec3& x) = 0;
	virtual Vec3 velAfterCollision(const Vec3& vel, const Vec3& n, float depth) = 0;
	virtual Vec3 getNormal(const Vec3& x) = 0;
};

struct Box : public Boundary {
	Vec3 c;
	float width; //x
	float height; //y
	float depth; //z
	Vec3 ext;

	Box(Vec3 c = Vec3{}, float width = 0.7f, float height = 0.7f, float depth = 0)
		: c{ c }, width{ width }, height{ height }, depth{ depth } {
		ext = Vec3{ width, height, depth };
	}

	Vec3 xToLocal(const Vec3& x);
	virtual Vec3 getContactPoint(const Vec3& x);
	virtual float F(const Vec3& x);
	virtual float getDepth(const Vec3& x);
	virtual Vec3 getSurfaceNormal(const Vec3& x);
	virtual Vec3 velAfterCollision(const Vec3& vel, const Vec3& n, float depth);
	virtual Vec3 getNormal(const Vec3& x);
};

struct Sphere : public Boundary {
	Vec3 c;
	float r;

	Sphere(Vec3 c = Vec3{}, float r = Const::borderR) : c{ c }, r{ r } {}

	virtual float F(const Vec3& x);
	virtual Vec3 getContactPoint(const Vec3& x);
	virtual float getDepth(const Vec3& x);
	virtual Vec3 getSurfaceNormal(const Vec3& x);
	virtual Vec3 velAfterCollision(const Vec3& vel, const Vec3& n, float depth);
	virtual Vec3 getNormal(const Vec3& x);
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

float B(float x, float y);

#endif // !SPHHELPER_HPP



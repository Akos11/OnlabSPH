#include "Sph.hpp"

void SPH_Simulator::densityAndPresureStep() {

	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles.particles[i];
		std::vector<Particle *> neighbors = particles.spatialQuery(p);

		p->rho = calcMassDensity(p, neighbors);
		p->p = calcPressure(p);
	}
}

float SPH_Simulator::calcMassDensity(const Particle * p, const std::vector<Particle *>& neighbors) {
	float rho = 0.0f;
	Vec3 r = p->pos;

	for (auto n : neighbors) {
		Vec3 rj = n->pos;
		rho += Const::particleM * W_def(r - rj, Const::h);
	}
	return rho;
}

float SPH_Simulator::calcPressure(const Particle * p) {
	float rho = p->rho;
	return Const::k * (rho - Const::rho0);
}

void SPH_Simulator::forcesStep() {

	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles.particles[i];
		std::vector<Particle *> neighbors = particles.spatialQuery(p);

		internalForcesStep(p, neighbors);
		externalForcesStep(p, neighbors);

		p->F = p->Fint + p->Fext;
	}
}

void SPH_Simulator::internalForcesStep(Particle * p, const std::vector<Particle *>& neighbors) {

	p->Fint = Vec3{};
	//lc pressure force
	p->Fint += calcPressureForce(p, neighbors);
	//calc viscosity force
	p->Fint += calcViscosityForce(p, neighbors);
}

Vec3 SPH_Simulator::calcPressureForce(const Particle * particle, const std::vector<Particle *>& neighbors) {
	float rho = particle->rho;
	float p = particle->p;
	Vec3 r = particle->pos;

	Vec3 Fpressure = Vec3{};
	for (auto n : neighbors) {
		if (n != particle) {
			Vec3 rj = n->pos;
			float rhoj = n->rho;
			float pj = n->p;

			Fpressure += ((p / (rho*rho)) + (pj / (rhoj*rhoj))) * Const::particleM * W_press_grad(r - rj, Const::h);
		}
	}

	return (-rho) * Fpressure;
}

Vec3 SPH_Simulator::calcViscosityForce(const Particle * particle, const std::vector<Particle *>& neighbors) {
	Vec3 vel = particle->getVelocityAtT();
	Vec3 r = particle->pos;

	Vec3 Fviscosity = Vec3{};
	for (auto n : neighbors) {
		if (n != particle) {
			Vec3 rj = n->pos;
			Vec3 velj = n->getVelocityAtT();
			float rhoj = n->rho;

			Fviscosity += (Const::particleM / rhoj) * W_visc_lapl(r - rj, Const::h) * (velj - vel);
		}
	}

	return Const::mu * Fviscosity;
}

void SPH_Simulator::externalForcesStep(Particle * p, const std::vector<Particle *>& neighbors) {
	Vec3 Fgravity = calcGravityForce(p, neighbors);
	Vec3 Fsurface = calcTensionForce(p, neighbors);

	p->Fext = Fgravity + Fsurface;
}

Vec3 SPH_Simulator::calcGravityForce(const Particle * p, const std::vector<Particle *>& neighbors) {
	return p->rho * Const::g;
}

Vec3 SPH_Simulator::calcTensionForce(const Particle * p, const std::vector<Particle *>& neighbors) {
	Vec3 n = calcInwardSurfaceNormal(p, neighbors);

	Vec3 Fsurface = Vec3{};
	if (n.len() >= Const::threshold){
		Fsurface = (-Const::surfTension) * calcSurfaceColorFieldLapl(p, neighbors) * (n / n.len());
	}

	return Fsurface;
}

Vec3 SPH_Simulator::calcInwardSurfaceNormal(const Particle * p, const std::vector<Particle *>& neighbors) {
	Vec3 r = p->pos;

	Vec3 normal = Vec3{};
	for (auto n : neighbors) {
		float rhoj = n->rho;
		Vec3 rj = n->pos;

		normal += (Const::particleM / rhoj) * W_def_grad(r - rj, Const::h);
	}

	return normal;
}

float SPH_Simulator::calcSurfaceColorFieldLapl(const Particle * p, const std::vector<Particle *>& neighbors) {
	Vec3 r = p->pos;

	float colorFieldLapl = 0.0f;
	for (auto n : neighbors) {
		float rhoj = n->rho;
		Vec3 rj = n->pos;

		colorFieldLapl += (Const::particleM / rhoj) * W_def_lapl(r - rj, Const::h);
	}

	return colorFieldLapl;
}

void SPH_Simulator::collisionHandlingStep() {
	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles.particles[i];

		//float boundaryBoxF = boundingBox.F(p->pos);
		float boundaryF = boundary->F(p->pos);

		if (boundaryF > 0.0f) {
			Vec3 contactPoint = boundary->getContactPoint(p->pos);
			p->pos = contactPoint;// -0.01f * (contactPoint / contactPoint.len());// -randFloatBtw(0.0f, 0.01) * (contactPoint / contactPoint.len());
			p->currVel = boundary->velAfterCollision(p->currVel, boundary->getSurfaceNormal(p->pos), boundary->getDepth(p->pos));

			//std::cout << contactPoint;
			//std::cout << boundary->getSurfaceNormal(p->pos) << std::endl;
		}

		p->currVel = p->getVelocityAtT();
	}
}

void SPH_Simulator::boundaryForceStep() {
	for (int i = 0; i < Const::particleNum; ++i) {
		Particle * p = particles.particles[i];
		std::vector<Particle *> borderNeighbors = particles.spatialQueryBorder(p);

		Vec3 Fak = Vec3{};
		for (auto k : borderNeighbors) {
			Vec3 xak = (k->pos) - (p->pos);
			//Vec3 surfNorm = k->pos / k->pos.len();
			Vec3 surfNorm = boundary->getNormal(Const::box ? p->pos : k->pos);

			float y = xak.dot(surfNorm);
			float x = sqrt(xak.len() * xak.len() - y*y);

			Fak += -(Const::particleM / (Const::particleM + Const::particleM)) * B(x, y) * surfNorm;
		}

		p->F += Fak;
	}
}

//PUBLIC

void SPH_Simulator::init() {
	particles = Particles{};
	//		boundary = new Sphere{};
	if (Const::box)
		boundary = new Box{};
	else
		boundary = new Sphere{};

	integrator = LeapFrogIntegrator{};

	particles.init();

	if (Const::box)
		particles.initBorderBox();
	else
		particles.initBorder();
}

void SPH_Simulator::simulationStep() {
	densityAndPresureStep();

	forcesStep();
	boundaryForceStep();


	integrator.stepIntegrator(particles.particles);

	collisionHandlingStep();

	particles.updateSpatialHashing();
}

void SPH_Simulator::simulate(float t) {
	while (integrator.t < t)
		simulationStep();
}

const std::vector<Particle *>& SPH_Simulator::getParticles() {
	return particles.particles;
}

const std::vector<Particle *>& SPH_Simulator::getBorderParticles() {
	return particles.borderParticles;
}

const std::vector<Particle *> SPH_Simulator::testBorderSpatialQuery(unsigned int i) {
	std::vector<Particle *> neighs = particles.spatialQueryBorder(particles.particles[i]);
	return neighs;
}

const std::vector<Particle *> SPH_Simulator::testSpatialQuery(unsigned int i) {
	std::vector<Particle *> neighs = particles.spatialQuery(particles.particles[i]);
	return neighs;
}

void SPH_Simulator::clear() {
	for (int i = 0; i < Const::particleNum; ++i) {
		delete particles.particles[i];
	}

	for (int i = 0; i < Const::borderParticleNum; ++i) {
		delete particles.borderParticles[i];
	}

	delete boundary;
}
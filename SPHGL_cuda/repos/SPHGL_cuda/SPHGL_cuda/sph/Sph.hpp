#ifndef SPH_HPP
#define SPH_HPP

#include "SphHelper.hpp"
#include "Const.hpp"
#include "Vec.hpp"

#include "..\kernel.h"

class SPH_Simulator {

	Particles particles;
	Boundary * boundary;
	LeapFrogIntegrator integrator;


	void densityAndPresureStep();
	float calcMassDensity(const Particle * p, const std::vector<Particle *>& neighbors);
	float calcPressure(const Particle * p);
	void forcesStep();
	void internalForcesStep(Particle * p, const std::vector<Particle *>& neighbors);
	Vec3 calcPressureForce(const Particle * particle, const std::vector<Particle *>& neighbors);
	Vec3 calcViscosityForce(const Particle * particle, const std::vector<Particle *>& neighbors);
	void externalForcesStep(Particle * p, const std::vector<Particle *>& neighbors);
	Vec3 calcGravityForce(const Particle * p, const std::vector<Particle *>& neighbors);
	Vec3 calcTensionForce(const Particle * p, const std::vector<Particle *>& neighbors);
	Vec3 calcInwardSurfaceNormal(const Particle * p, const std::vector<Particle *>& neighbors);
	float calcSurfaceColorFieldLapl(const Particle * p, const std::vector<Particle *>& neighbors);
	void collisionHandlingStep();
	void boundaryForceStep();


public:
	SPH_Simulator() {}
	void init();
	void simulationStep();
	void simulate(float t);
	const std::vector<Particle *>& getParticles();
	const std::vector<Particle *>& getBorderParticles();
	const std::vector<Particle *> testBorderSpatialQuery(unsigned int i);
	const std::vector<Particle *> testSpatialQuery(unsigned int i);
	void clear();

};

////END OF Simulation

#endif // !SPHHELPER_HPP
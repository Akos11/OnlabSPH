// SphereView.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "stdafx.h"
#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "SphHelper.hpp"
#include "Const.hpp"
#include "Vec.hpp"

int width = 512;
int height = 512;

bool keysPressed[256];
float radius = 0.02f;

////Simulation

class SPH_Simulator {

	Particles particles;
//	Box sphere;
	Sphere sphere;
	LeapFrogIntegrator integrator;


	void densityAndPresureStep() {

		for (int i = 0; i < Const::particleNum; ++i) {
			Particle * p = particles.particles[i];
			std::vector<Particle *> neighbors = particles.spatialQuery(p);

			//calc mass-density
			p->rho = calcMassDensity(p, neighbors);
			//calc pressure
			p->p = calcPressure(p);
		}
	}

	float calcMassDensity(const Particle * p, const std::vector<Particle *>& neighbors) {
		float rho = 0.0f;
		Vec3 r = p->pos;

		//		std::cout << neighbors.size() << std::endl;

		for (auto n : neighbors) {
			Vec3 rj = n->pos;
			rho += Const::particleM * W_def(r - rj, Const::h);
		}
		return rho;
	}

	float calcPressure(const Particle * p) {
		float rho = p->rho;
		return Const::k * (rho - Const::rho0);
	}

	void forcesStep() {

		for (int i = 0; i < Const::particleNum; ++i) {
			Particle * p = particles.particles[i];
			std::vector<Particle *> neighbors = particles.spatialQuery(p);

			internalForcesStep(p, neighbors);
			externalForcesStep(p, neighbors);
			
			p->F = p->Fint + p->Fext;
		}
	}

	void internalForcesStep(Particle * p, const std::vector<Particle *>& neighbors) {

		p->Fint = Vec3{};
		//calc pressure force
		p->Fint += calcPressureForce(p, neighbors);
		//calc viscosity force
		p->Fint += calcViscosityForce(p, neighbors);
	}

	Vec3 calcPressureForce(const Particle * particle, const std::vector<Particle *>& neighbors) {
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

	Vec3 calcViscosityForce(const Particle * particle, const std::vector<Particle *>& neighbors) {
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

	void externalForcesStep(Particle * p, const std::vector<Particle *>& neighbors) {
		Vec3 Fgravity = calcGravityForce(p, neighbors);
		Vec3 Fsurface = calcTensionForce(p, neighbors);

		p->Fext = Fgravity + Fsurface;
	}

	Vec3 calcGravityForce(const Particle * p, const std::vector<Particle *>& neighbors) {
		return p->rho * Const::g;
	}

	Vec3 calcTensionForce(const Particle * p, const std::vector<Particle *>& neighbors) {
		Vec3 n = calcInwardSurfaceNormal(p, neighbors);

		Vec3 Fsurface = Vec3{};
		if (n.len() >= Const::threshold){
			Fsurface = (-Const::surfTension) * calcSurfaceColorFieldLapl(p, neighbors) * (n / n.len());
		}

		return Fsurface;
	}

	Vec3 calcInwardSurfaceNormal(const Particle * p, const std::vector<Particle *>& neighbors) {
		Vec3 r = p->pos;

		Vec3 normal = Vec3{};
		for (auto n : neighbors) {
			float rhoj = n->rho;
			Vec3 rj = n->pos;

			normal += (Const::particleM / rhoj) * W_def_grad(r - rj, Const::h);
		}

		return normal;
	}

	float calcSurfaceColorFieldLapl(const Particle * p, const std::vector<Particle *>& neighbors) {
		Vec3 r = p->pos;

		float colorFieldLapl = 0.0f;
		for (auto n : neighbors) {
			float rhoj = n->rho;
			Vec3 rj = n->pos;

			colorFieldLapl += (Const::particleM / rhoj) * W_def_lapl(r - rj, Const::h);
		}

		return colorFieldLapl;
	}

	void collisionHandlingStep() {
		for (int i = 0; i < Const::particleNum; ++i) {
			Particle * p = particles.particles[i];

			//float boundaryBoxF = boundingBox.F(p->pos);
			float sphereF = sphere.F(p->pos);

			if (sphereF > 0.0f) {
				Vec3 contactPoint = sphere.getContactPoint(p->pos);
				p->pos = contactPoint -randFloatBtw(0.0f, 0.01) * (contactPoint / contactPoint.len());
				p->currVel = sphere.velAfterCollision(p->currVel, sphere.getSurfaceNormal(p->pos), sphere.getDepth(p->pos));
			}

			p->currVel = p->getVelocityAtT();
		}
	}
	

public:
	SPH_Simulator() {}

	void init() {
		particles = Particles{};
		//boundingBox = Box{};
		sphere = Sphere{};
//		sphere = Box{};
		integrator = LeapFrogIntegrator{};

		particles.init();
	}

	void simulationStep() {
		densityAndPresureStep();

		forcesStep();


		integrator.stepIntegrator(particles.particles);
		
		collisionHandlingStep();

		particles.updateSpatialHashing();
	}

	void simulate(float t) {
		while (integrator.t < t)
			simulationStep();
	}

	const std::vector<Particle *>& getParticles() {
		return particles.particles;
	}


};

////END OF Simulation

SPH_Simulator sim = SPH_Simulator{};

void initOpenGL() {
	glewExperimental = GL_TRUE;
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		std::cerr << "Error: " << glewGetErrorString(err) << std::endl;
	}
	else {
		if (GLEW_VERSION_3_0)
		{
			std::cout << "Driver supports OpenGL 3.0\nDetails:" << std::endl;
			std::cout << "  Using GLEW " << glewGetString(GLEW_VERSION) << std::endl;
			std::cout << "  Vendor: " << glGetString(GL_VENDOR) << std::endl;
			std::cout << "  Renderer: " << glGetString(GL_RENDERER) << std::endl;
			std::cout << "  Version: " << glGetString(GL_VERSION) << std::endl;
			std::cout << "  GLSL: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
		}
	}

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void DrawCircle(float cx, float cy, float r, int num_segments)
{
	glBegin(GL_TRIANGLE_FAN);
	for (int ii = 0; ii < num_segments; ii++)
	{
		float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

		float x = r * cosf(theta);//calculate the x component
		float y = r * sinf(theta);//calculate the y component

		glVertex2f(x + cx, y + cy);//output vertex

	}
	glEnd();
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);

//	simulationStep();

	glColor3f(0.17f, 0.4f, 0.6f);
	//visualizationStep();
//	for (int i = 0; i < particleNumber; i++)
//	{
//		DrawCircle(particles[i].v[0],particles[i].v[1],radius,30);
//	}

	sim.simulationStep();

	const std::vector<Particle *> prtcls = sim.getParticles();
	for (auto p : prtcls) {
		Vec3 r = p->pos;

		DrawCircle(r.x, r.y, radius, 30);
	}



	glEnable(GL_DEPTH_TEST);
	glutSwapBuffers();
}

void idle() {
	glutPostRedisplay();
}

void keyDown(unsigned char key, int x, int y) {
	keysPressed[key] = true;
}

void keyUp(unsigned char key, int x, int y) {
	keysPressed[key] = false;
	switch (key) {
	case 'r':
//		resetSimulation();
		sim.init();
		break;
	}
/*
	case 'd':
		visualizationMethod = 0;
		break;
	case 'v':
		visualizationMethod = 1;
		break;
	case 'p':
		visualizationMethod = 2;
		break;

	case '1':
		densityColor.s[0] = densityColor.s[1] = densityColor.s[2] = densityColor.s[3] = 1.0f;
		break;

	case '2':
		densityColor.s[0] = 1.0f;
		densityColor.s[1] = densityColor.s[2] = densityColor.s[3] = 0.0f;
		break;

	case '3':
		densityColor.s[1] = 1.0f;
		densityColor.s[0] = densityColor.s[2] = densityColor.s[3] = 0.0f;
		break;

	case '4':
		densityColor.s[2] = 1.0f;
		densityColor.s[0] = densityColor.s[1] = densityColor.s[3] = 0.0f;
		break;

	case 27:
		exit(0);
		break;
	}*/
}

int mX, mY;

void mouseClick(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON)
		if (state == GLUT_DOWN) {
			mX = x;
			mY = y;
		}
}

void mouseMove(int x, int y) {
	/*force.s[0] = (float)(x - mX) * 10;
	force.s[1] = -(float)(y - mY) * 10;
	addForce(mX, height - mY, force);
	//addForce(256, 256, force);*/
	mX = x;
	mY = y;
}

void reshape(int newWidth, int newHeight) {
	width = newWidth;
	height = newHeight;
	glViewport(0, 0, width, height);
}


int main(int argc, char* argv[]) {
	glutInit(&argc, argv);
	glutInitContextVersion(3, 0);
	glutInitContextFlags(GLUT_CORE_PROFILE | GLUT_DEBUG);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(width, height);
	glutCreateWindow("SphereView");

	initOpenGL();

	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyDown);
	glutKeyboardUpFunc(keyUp);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);

	// OpenCL processing
//	initSimulation();
	
//	partics.init();

	sim.init();
	//sim.simulate(1);

	glutMainLoop();
	return(0);
}


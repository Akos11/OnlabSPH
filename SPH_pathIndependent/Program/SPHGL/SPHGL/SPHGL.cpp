// SphereView.cpp : This file contains the 'main' function. Program execution begins and ends here.
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
float radius = 0.015f;

int mode = 0;

////Simulation

class SPH_Simulator {

	Particles particles;
//	Box sphere;
//	Sphere boundary;
	Boundary * boundary;
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

	void boundaryForceStep() {
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
	

public:
	SPH_Simulator() {}

	void init() {
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

	void simulationStep() {
		densityAndPresureStep();

		forcesStep();
		boundaryForceStep();


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

	const std::vector<Particle *>& getBorderParticles() {
		return particles.borderParticles;
	}

	const std::vector<Particle *> testBorderSpatialQuery(unsigned int i) {
		std::vector<Particle *> neighs =  particles.spatialQueryBorder(particles.particles[i]);
		return neighs;
	}

	void clear() {
		for (int i = 0; i < Const::particleNum; ++i) {
			delete particles.particles[i];
		}

		for (int i = 0; i < Const::borderParticleNum; ++i) {
			delete particles.borderParticles[i];
		}

		delete boundary;
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

	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	//gluLookAt(1.0, 1.0, -1.0, 0, 0, 0, 1, 0, 0);
}

void DrawCircle(float cx, float cy, float r, int num_segments, Vec3* velocity = NULL, Vec3* color = NULL)
{
	glBegin(GL_TRIANGLE_FAN);
	float length = -1;
	if (velocity != NULL) {
		length = velocity->len();
		
	}
	if (color != NULL) {
		glColor3f(color->x, color->y, color->z);
	}
	else if (length != -1){
		glColor3f(length, 0.1f, 0.1f);
	}
	for (int ii = 0; ii < num_segments; ii++)
	{
		float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

		float x = r * cosf(theta);//calculate the x component
		float y = r * sinf(theta);//calculate the y component

		if (velocity != NULL)
		{
			Vec3 temp = Vec3(x, y, 0);
			temp = length * (temp / r);
			if (fabs(temp.x - velocity->x) < 0.1 && fabs(temp.y - velocity->y) < 0.1)
			{
				float tempf = 12.0f * length;
				if (tempf < 1) tempf = 1;
				x *= tempf;
				y *= tempf;
			}
		}
		glVertex2f(x + cx, y + cy);//output vertex

	}
	glEnd();
}
///uj3d
double eyeX = 1.0;
double eyeY = 1.0;
double eyeZ = 1.0;
///uj3d

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	///uj3d
	glLoadIdentity();
	gluPerspective(45.0f, 1.0, 0.01, 200.0);
	///uj3d

	//glViewport(0, 0, 1000, 1000);
//	simulationStep();

	//glColor3f(0.5f, 0.5f, 0.5f);

	//DrawCircle(0, 0, 0.8f, 60);

	glColor3f(0.17f, 0.4f, 0.6f);
	//visualizationStep();
//	for (int i = 0; i < particleNumber; i++)
//	{
//		DrawCircle(particles[i].v[0],particles[i].v[1],radius,30);
//	}

	sim.simulationStep();
	///uj3d
	gluLookAt(eyeX, eyeY, eyeZ, 0, 0, 0.0, 0.0, 1.0, 0.0);
	///uj3d

	if (!Const::DDD) {
		glColor3f(0.0f, 1.0f, 0.0f);
		const std::vector<Particle *> border = sim.getBorderParticles();
		for (auto p : border) {
			Vec3 r = p->pos;
			//DrawCircle(r.x, r.y, radius, 7);
			///uj3d
			glTranslatef(r.x, r.y,0);
			glutWireSphere(radius, 20, 10);
			glTranslatef(-r.x, -r.y, 0);
			///uj3d
		}
	}
	glColor3f(0.0f, 0.0f, 1.0f);

	const std::vector<Particle *> prtcls = sim.getParticles();
	for (auto p : prtcls) {
		Vec3 r = p->pos;
		float v = p->currVel.len();
		//glColor3f(v, v, 1.0f);
	/*	if (mode == 1)
			DrawCircle(r.x, r.y, radius, 7,&(p->currVel), &(p->color));
		else if (mode == 2)
			DrawCircle(r.x, r.y, radius, 7, &(p->currVel));
		else
			DrawCircle(r.x, r.y, radius, 7);*/
		///uj3d
		glTranslatef(r.x, r.y, r.z);
		glutWireSphere(radius, 20, 10);
		glTranslatef(-r.x, -r.y, -r.z);
		///uj3d

	}
	
	

	/*
	
	glColor3f(1.0f, 0.0f, 0.0f);
	std::vector<Particle *> testBorderNeighbors = sim.testBorderSpatialQuery(55);
	for (auto p : testBorderNeighbors) {
		Vec3 r = p->pos;
		DrawCircle(r.x, r.y, radius, 10);
	}
	*/	
	


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
		sim.clear();
		sim.init();
		break;
	
	case 'v':
		mode = (mode+1) % 3;
		break;
	case 'w':
		eyeZ += 0.05;
		break;
	case 's':
		eyeZ -= 0.05;
		break;
	case 'd':
		eyeX += 0.05;
		break;
	case 'a':
		eyeX -= 0.05;
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


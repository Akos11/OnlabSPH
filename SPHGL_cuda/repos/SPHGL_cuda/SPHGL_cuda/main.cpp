// SphereView.cpp : This file contains the 'main' function. Program execution begins and ends here.
//

//#include "stdafx.h"
#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "sph\SphHelper.hpp"
#include "sph\Const.hpp"
#include "sph\Vec.hpp"
#include "sph\Sph.hpp"

#include "kernel.h"

int width = 512;
int height = 512;

bool keysPressed[256];
float radius = 0.015f;

int mode = 0;

SPH_Simulator sim = SPH_Simulator{};

std::vector<Vec3> cudaParticles = std::vector<Vec3>();

int nFrames = 0;

void cudaUpdateParticles() {
	float3* positionBufferCPU = new float3[Const::particleNum];

	cuda::retrievePositionData(positionBufferCPU);
	//cuda::testSpatialQuery(positionBufferCPU, 342);
	
	cudaParticles.clear();
	for (int i = 0; i < Const::particleNum; ++i) {
		float3 pos = positionBufferCPU[i];
		
		/*if (i == 134)
			printf("%f %f\n", pos.x, pos.y);*/
		
		cudaParticles.push_back(Vec3{ pos.x, pos.y, pos.z });
	}

	delete[] positionBufferCPU;
}


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
		float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);	//get the current angle

		float x = r * cosf(theta);	//calculate the x component
		float y = r * sinf(theta);	//calculate the y component

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
		glVertex2f(x + cx, y + cy);	//output vertex

	}
	glEnd();
}

void displayParticles() {
	const std::vector<Particle *> prtcls = sim.getParticles();
	for (auto p : prtcls) {
		Vec3 r = p->pos;
		float v = p->currVel.len();
		if (mode == 1)
			DrawCircle(r.x, r.y, radius, 7, &(p->currVel), &(p->color));
		else if (mode == 2)
			DrawCircle(r.x, r.y, radius, 7, &(p->currVel));
		else
			DrawCircle(r.x, r.y, radius, 7);
	}

	if (!Const::DDD) {
		glColor3f(0.0f, 1.0f, 0.0f);
		const std::vector<Particle *> border = sim.getBorderParticles();
		for (auto p : border) {
			Vec3 r = p->pos;
			DrawCircle(r.x, r.y, radius, 7);
		}
	}
}

void cudaDisplayParticles() {
	cudaUpdateParticles();

	//cuda::simulationStep();

	for (auto p : cudaParticles) {
		Vec3 r = p;

		//std::cout << r.z;
		
		//DrawCircle(r.x, r.y, radius, 7);
		///uj3d
		glTranslatef(r.x, r.y, r.z);
		glutWireSphere(radius, 20, 10);
		glTranslatef(-r.x, -r.y, -r.z);
	}

	glColor3f(0.0f, 1.0f, 0.0f);
}

double eyeX = 1.0;
double eyeY = 1.0;
double eyeZ = 1.0;

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);

	glLoadIdentity();
	gluPerspective(45.0f, 1.0, 0.01, 200.0);

	glColor3f(0.17f, 0.4f, 0.6f);

	//sim.simulationStep();

	glColor3f(0.0f, 0.0f, 1.0f);
	
	//Display Particles

	//displayParticles();
	gluLookAt(eyeX, eyeY, eyeZ, 0, 0, 0.0, 0.0, 1.0, 0.0);
	cudaDisplayParticles();

	/*const std::vector<Particle *> neigh = sim.testSpatialQuery(322);

	for (auto p : neigh) {
		Vec3 r = p->pos;

		DrawCircle(r.x, r.y, radius, 7);
	}*/

	glColor3f(0.0f, 1.0f, 0.0f);
	/*const std::vector<Particle *> border = sim.getBorderParticles();
	for (auto p : border) {
		Vec3 r = p->pos;
		glTranslatef(r.x, r.y, 0);
		glutWireSphere(radius, 20, 10);
		glTranslatef(-r.x, -r.y, 0);
	}*/
	glutWireSphere(0.8, 20, 10);
	glEnable(GL_DEPTH_TEST);
	glutSwapBuffers();
}
float speed = 1.0;
void idle() {
	nFrames += 1;

	if (nFrames % 1000 == 0)
		printf("%d\n", nFrames);

	long time = glutGet(GLUT_ELAPSED_TIME); // elapsed time since the start of the program

	static float tend = 0;
	const float dt = 0.01f;
	float tstart = tend;
	tend = time / 1000.0f;
	//glutPostRedisplay();					// redraw the scene

	for (float t = tstart; t < tend; t += dt) {
		float Dt = fmin(dt, tend - t);
		
		//cuda::simulationStep(Dt);
		//glutPostRedisplay();
	}
	cuda::simulationStep(dt/speed);
	glutPostRedisplay();					// redraw the scene
}

void keyDown(unsigned char key, int x, int y) {
	keysPressed[key] = true;
}

void keyUp(unsigned char key, int x, int y) {
	keysPressed[key] = false;
	switch (key) {
	case 'r':
		sim.clear();
		sim.init();
		break;

	case 'v':
		mode = (mode + 1) % 3;
		break;

	case 'k':
		speed -= 0.5f;
		break;
	case 'j':
		speed += 0.5f;
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
	
	sim.init();

	//CUDA
	cuda::initSimulation();

	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyDown);
	glutKeyboardUpFunc(keyUp);
	glutMouseFunc(mouseClick);
	glutMotionFunc(mouseMove);


	

	glutMainLoop();
	return(0);
}


//int main(int argc, char* argv[]) {
//	//CUDA
//	cuda::initSimulation();
//	
//	int i = 0;
//	printf("started\n");
//	while (i++ < 10000) {
//		cuda::simulationStep(0.01f);
//
//		if (i % 1000 == 0)
//			printf("%d\n", i);
//	}
//
//	return(0);
//}
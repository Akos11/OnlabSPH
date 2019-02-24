// SphereView.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>

#define PI 3.14159
int width = 512;
int height = 512;

bool keysPressed[256];
int particleNumber = 2000;
float radius = 0.03f;



struct vec3 {
	float v[3];


	vec3(float x = 0, float y = 0, float z = 0) { v[0] = x; v[1] = y; v[2] = z; }	
	vec3 operator*(float s) {
		vec3 res(v[0] * s, v[1] * s, v[2] * s);
		return res;
	}
};

vec3 * particles = new vec3[particleNumber];

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
void initSimulation() {
	for (size_t i = 0; i < particleNumber; i++)
	{
		//std::cout << (float)rand() / RAND_MAX;
		particles[i] = vec3(((float)rand() / RAND_MAX * 2) - 1, ((float)rand() / RAND_MAX * 2) - 1, 0.0f);
	}

}
void resetSimulation() {

	initSimulation();
}
void simulationStep() {
	for (size_t i = 0; i < particleNumber; i++)
	{
		if (particles[i].v[1] > -0.95f)
		particles[i].v[1] -= 0.001f;
	}
}
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	simulationStep();
	glColor3f(0.17f, 0.4f, 0.6f);
	//visualizationStep();
	for (int i = 0; i < particleNumber; i++)
	{
		DrawCircle(particles[i].v[0],particles[i].v[1],radius,30);
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
		resetSimulation();
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
	initSimulation();

	glutMainLoop();
	return(0);
}


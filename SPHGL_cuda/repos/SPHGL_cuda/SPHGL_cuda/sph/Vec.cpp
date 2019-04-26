//#include "stdafx.h"
#include <cmath>
#include "Vec.hpp"

Vec2 Vec2::operator+(const Vec2& other) const {
	return Vec2{ this->x + other.x, this->y + other.y };
}

Vec2& Vec2::operator+=(const Vec2& other) {
	this->x += other.x;
	this->y += other.y;

	return *this;
}

Vec2 Vec2::operator-(const Vec2& other) const {
	return Vec2{ this->x - other.x, this->y - other.y };
}

Vec2& Vec2::operator-=(const Vec2& other) {
	this->x -= other.x;
	this->y -= other.y;

	return *this;
}

float Vec2::dot(const Vec2& other) const {
	return this->x * other.x + this->y * other.y;
}

std::ostream& operator<<(std::ostream& os, const Vec2& v) {
	os << "(" << v.x << ", " << v.y << ")" << std::endl;
	return os;
}

Vec3 Vec3::operator+(const Vec3& other) const {
	return Vec3{ this->x + other.x, this->y + other.y, this->z + other.z };
}

Vec3& Vec3::operator+=(const Vec3& other) {
	this->x += other.x;
	this->y += other.y;
	this->z += other.z;

	return *this;
}

Vec3 Vec3::operator-(const Vec3& other) const {
	return Vec3{ this->x - other.x, this->y - other.y, this->z - other.z };
}

Vec3& Vec3::operator-=(const Vec3& other) {
	this->x -= other.x;
	this->y -= other.y;
	this->z -= other.z;

	return *this;
}

Vec3 Vec3::operator/(const float f) const {
	return Vec3{ this->x / f, this->y / f, this->z / f };
}

float Vec3::dot(const Vec3& other) const {
	return this->x * other.x + this->y * other.y + this->z * other.z;
}

float Vec3::len() const {
	return sqrt(x*x + y*y + z*z);
}
float Vec3::dist(const Vec3& other) const {
	return (*this - other).len();
}

Vec3 Vec3::operator-() const {
	return Vec3{ -this->x, -this->y, -this->z };
}

Vec3 Vec3::callFuncForComponents(std::function<float(float)> f) {
	return Vec3{ f(this->x), f(this->y), f(this->z) };
}

float Vec3::getMaxComponent() {
	//float yOrZ = this->y > this->z ? this->y : this->z;
	//return this->x > yOrZ ? this->x : yOrZ;
	return this->x > this->y ? this->x : this->y;
}
Vec3 callFuncForComponents(std::function<float(float, float)> f, const Vec3& v1, const Vec3& v2) {
	return Vec3{ f(v1.x, v2.x), f(v1.y, v2.y), f(v1.z, v2.z) };
}

std::ostream& operator<<(std::ostream& os, const Vec3& v) {
	os << v.x << " " << v.y << " " << std::endl;
	return os;
}

Vec3 operator*(float f, const Vec3& v) {
	return Vec3{ f*v.x, f*v.y, f*v.z };
}


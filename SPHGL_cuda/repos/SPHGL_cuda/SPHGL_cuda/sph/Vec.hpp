#ifndef VEC_HPP
#define VEC_HPP

#include <iostream>
#include <functional>

struct Vec2 {

	float x;
	float y;

	Vec2(float x = 0.0, float y = 0.0) : x{ x }, y{ y } {}

	Vec2 operator+(const Vec2& other) const;
	Vec2& operator+=(const Vec2& other);
	Vec2 operator-(const Vec2& other) const;
	Vec2& operator-=(const Vec2& other);

	float dot(const Vec2& other) const;

};


std::ostream& operator<<(std::ostream& os, const Vec2& v);

struct Vec3 {

	float x;
	float y;
	float z;

	Vec3(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x{ x }, y{ y }, z{ z } {}

	Vec3 operator+(const Vec3& other) const;
	Vec3& operator+=(const Vec3& other);
	Vec3 operator-(const Vec3& other) const;
	Vec3& operator-=(const Vec3& other);

	Vec3 operator/(const float f) const;

	float dot(const Vec3& other) const;

	float dist(const Vec3& other) const;
	float len() const;

	Vec3 operator-() const;

	Vec3 callFuncForComponents(std::function<float(float)> f);

	float getMaxComponent();
};


std::ostream& operator<<(std::ostream& os, const Vec3& v);

Vec3 callFuncForComponents(std::function<float(float, float)> f, const Vec3& v1, const Vec3& v2);

Vec3 operator*(float f, const Vec3& v);

#endif


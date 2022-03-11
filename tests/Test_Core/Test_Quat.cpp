#include "gtest/gtest.h"
#include "Quat.h"
#include "Matrix.h"
#include <cstdlib>
#include <ctime>

#define cos_angle(angle) cos(double(angle * 3.1415926f / 180.0))
#define sin_angle(angle) sin(double(angle * 3.1415926f / 180.0))
#define radian(x) x * 3.1415926f / 180.0
using namespace dyno;

Mat4f getTransform(Vec3f & T, Vec3f& R, Vec3f& S)
{
	// TODO:改为四元数
	Mat4f translation = Mat4f(
		1, 0, 0, T[0],
		0, 1, 0, T[1],
		0, 0, 1, T[2],
		0, 0, 0, 1);
	
	// R[X,Y,Z] -> [Z,X,Y]轴

	double X = R[0];
	Mat4f rotation_x = Mat4f(
		1, 0, 0, 0,
		0, cos_angle(X), -sin_angle(X), 0,
		0, sin_angle(X), cos_angle(X), 0,
		0, 0, 0, 1);

	double Y = R[1];
	Mat4f rotation_y = Mat4f(
		cos_angle(Y), 0, sin_angle(Y), 0,
		0, 1, 0, 0,
		-sin_angle(Y), 0, cos_angle(Y), 0,
		0, 0, 0, 1);

	double Z = R[2];
	Mat4f rotation_z = Mat4f(
		cos_angle(Z), -sin_angle(Z), 0, 0,
		sin_angle(Z), cos_angle(Z), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	Mat4f scaling= Mat4f(
		S[0], 0, 0, 0,
		0, S[1], 0, 0,
		0, 0, S[2], 0,
		0, 0, 0, 1);

	return  translation * scaling * rotation_x * rotation_y * rotation_z;
}

void getQuat(Vec3f &T, Vec3f &R, float &S, Quat<float> &GlT, float &GlS, Quat<float> &GlR)
{

	// Scaling
	// Quat<Real> s(S[0], S[1], S[2], 0.f);
	// Translation
	Quat<Real> t(T[0], T[1], T[2], 0.f);
	// Rotate
	Quat<Real> q_x(radian(R[0]), Vec3f(1.f, 0.f, 0.f));
	Quat<Real> q_y(radian(R[1]), Vec3f(0.f, 1.f, 0.f));
	Quat<Real> q_z(radian(R[2]), Vec3f(0.f, 0.f, 1.f));
	
	GlT = GlT + GlS * GlR * t * GlR.conjugate();
	GlS = GlS * S;
	GlR = GlR * (q_x * q_y * q_z);

	// Quat<Real> tmp = GlR;
	// TODO:
	// GlR.normalize(); 
	// EXPECT_EQ(GlR == tmp, true);
}

Vec3f getVec3fM(Vec3f x, Mat4f g)
{
	Vec4f tmp_p(x[0], x[1], x[2], 1);
	
	tmp_p = g * tmp_p;
	// PT_f("tmp3", tmp_p[3], pId);
	x[0] = tmp_p[0] / tmp_p[3];
	x[1] = tmp_p[1] / tmp_p[3];
	x[2] = tmp_p[2] / tmp_p[3];

	return x;
}

Vec3f getVec3fQ(Vec3f x,  Quat<float> &GlT, float &GlS, Quat<float> &GlR)
{
	Quat<float> tmp_p(x[0], x[1], x[2], 0);
	
	tmp_p = GlS * GlR * tmp_p *  GlR.conjugate() + GlT;

	return Vec3f(tmp_p.x, tmp_p.y, tmp_p.z);
}

Vec3f getVec3fQInv(Vec3f x,  Quat<float> &GlT, float &GlS, Quat<float> &GlR)
{
	Quat<float> tmp_p(x[0], x[1], x[2], 0);
	
	tmp_p = (1. / GlS) * GlR.conjugate() * (tmp_p - GlT) *  GlR;

	return Vec3f(tmp_p.x, tmp_p.y, tmp_p.z);
}

const int N = 3;

TEST(Quat, func)
{
	//TEST Animation

	srand(time(0));
	Vec3f Translation[N] = {Vec3f(1.0,2.0,3.0), Vec3f(3.0,2.0,1.0), Vec3f(1.0,2.0,3.0)};
	// Vec3f Translation[N] = {Vec3f(0.0,0.0,0.0), Vec3f(0.0,0.0,0.0), Vec3f(0.0,0.0,0.0)};
	Vec3f Rotation[N] = {Vec3f(0,10,5), Vec3f(20,0,0), Vec3f(0,10,10)};
    // Vec3f Scaling[N] = {Vec3f(1.0,1.0,1.0), Vec3f(1.0,1.0,1.0), Vec3f(1.0,1.0,1.0)};
	float Scaling[N] = {4.0, 2.0, 3.0};
	

	Mat4f Transform[N];

	Quat<float> GlT[N];
	Quat<float> GlR[N];
	// Quat<float> GlS[N];
	float GlS[N];

	Transform[0] = getTransform(Translation[0], Rotation[0], Vec3f(Scaling[0]));
	for(int i = 1; i < N; ++i)
		Transform[i] = Transform[i-1] * getTransform(Translation[i], Rotation[i], Vec3f(Scaling[i]));

	GlT[0] = Quat<float>(0.f, 0.f, 0.f, 0.f);
    // GlS[0] = Quat<float>(1.f, 1.f, 1.f, 0.f);
	GlS[0] = 1.0f;
    GlR[0] = Quat<float>(0.f, 0.f, 0.f, 1.f);
	getQuat(Translation[0], Rotation[0], Scaling[0], GlT[0], GlS[0], GlR[0]);
	for(int i = 1; i < N; ++i)
	{
		GlT[i] = GlT[i - 1];
		GlS[i] = GlS[i - 1];
		GlR[i] = GlR[i - 1];
		getQuat(Translation[i], Rotation[i], Scaling[i], GlT[i], GlS[i], GlR[i]);
	}

	for(int i = 0; i < N; ++i)
	{
		Vec3f x(0.f, 0.f, 0.f);
		Vec3f xq = getVec3fQ(x, GlT[i], GlS[i], GlR[i]);
		Vec3f xm = getVec3fM(x, Transform[i]);

		// inverse
		Mat4f invM = Transform[i].inverse();
		Vec3f ym = getVec3fM(xm, invM);
		Vec3f yq = getVec3fQInv(xq, GlT[i], GlS[i], GlR[i]);		
		EXPECT_EQ( (xq - xm).norm() < 1e-5, true);
		EXPECT_EQ( (x - ym).norm() < 1e-5, true);
		EXPECT_EQ( (x - yq).norm() < 1e-5, true);		
		EXPECT_EQ( (x - yq).norm() < (x - ym).norm(), true);	
	}

	//Random
	for(int i = 0; i < N; ++i)
	{
		int count = 0;
		for(int _ = 0; _ < 10000; ++_)
		{
			Vec3f x;
			for (int j = 0; j < 3; ++j)
				x[j] = (rand() % 1000) / 10000000.0f;
			Quat<float> tmpR(Transform[i]);
			Vec3f xq = getVec3fQ(x, GlT[i], GlS[i], GlR[i]);
			// Vec3f xqm = getVec3fQ(x, GlT[i], GlS[i], tmpR);
			// Vec3f xmq = getVec3fM(x, GlR[i].toMatrix4x4());
			Vec3f xm = getVec3fM(x, Transform[i]);

			// inverse
			Mat4f invM = Transform[i].inverse();
			Vec3f ym = getVec3fM(xm, invM);
			Vec3f yq = getVec3fQInv(xq, GlT[i], GlS[i], GlR[i]);
			float disym = (x - ym).norm();
			float disyq = (x - yq).norm();
			count += (disyq < disym);
			EXPECT_EQ( (xq - xm).norm() < 1e-5, true);
			EXPECT_EQ( (x - ym).norm() < 1e-5, true);
			EXPECT_EQ( (x - yq).norm() < 1e-5, true);
			//EXPECT_EQ( disq < dism, true);	
		}
		printf("%f\%\n", count / 100.);
	}
}
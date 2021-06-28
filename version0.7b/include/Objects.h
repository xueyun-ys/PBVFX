#pragma once
#ifndef ____PBE_OBJECTS_H____
#define ____PBE_OBJECTS_H____

#include "Volume.h"
#include <cmath>

float PI = 3.1415926;

using namespace std;

namespace lux
{
	class ConstantField : public Volume<float>
	{

	public:

		ConstantField(float f) : F(f) {}
		~ConstantField() {}//;

		const float eval(const Vector& P) const { float base=1.0; /*std::cout<<"hello0!"<<std::endl;*/ return F; }
		const Vector grad(const Vector& P) const { return Vector(0, 0, 0); }


	private:

		double F;
	};

	class ObjectSphere : public Volume<float>
	{

	public:

		ObjectSphere(float r) : R(r) {}
		~ObjectSphere() {};

		const float eval(const Vector& X) const
		{
			//std::cout<<"Sphere!"<< R- X.magnitude()<<std::endl;
			return R - X.magnitude();
		}
		//Vector grad(const Vector& P) const { return Vector(0, 0, 0); }

	private:

		double R;
	};

	class ObjectEllipse : public Volume<float>
	{

	public:

		ObjectEllipse(Vector n, float r1, float r2) : N(n), R1(r1), R2(r2){}
		~ObjectEllipse() {};

		const float eval(const Vector& X) const
		{
			//std::cout<<"hello!"<<std::endl;
			float Z = X * N;
			float x_per = (X - Z * N).magnitude();
			x_per *= x_per;
			Z *= Z;
			float final = 1.0 - Z / (R1 * R1) - x_per / (R2 * R2);
			return final;
		}

	private:
		Vector N = Vector(1, 0, 0);
		float R1, R2;//R1:major; R2:minor
	};

	class ObjectTorus : public Volume<float>
	{

	public:

		ObjectTorus(Vector n, float r1, float r2) : N(n), R1(r1), R2(r2) {}
		~ObjectTorus() {};

		const float eval(const Vector& X)const
		{
			float x_per = (X - (X * N)*N).magnitude();
			x_per *= x_per;
			float final = 4.0*R1*R1*x_per-(     X.magnitude()* X.magnitude() + (R1 * R1) - (R2 * R2)        ) * (     X.magnitude() * X.magnitude() + (R1 * R1) - (R2 * R2)    );
			return final;
		}

	private:
		Vector N = Vector(1, 0, 0);
		float R1, R2;//R1:major; R2:minor
	};

	class ObjectBox : public Volume<float>
	{

	public:

		ObjectBox(float r, float q) : R(r), Q(q) {}
		~ObjectBox() {};

		const float eval(const Vector& X)const
		{

			float final = pow(R, Q*2) - pow((float)X.X(), Q*2) - pow((float)X.Y(), Q * 2) - pow((float)X.Z(), Q * 2);
			//float final = pow(R, Q*2) - pow((float)X[0], Q*2) - pow((float)X[1], Q * 2) - pow((float)X[2], Q * 2);
			return final;
		}

	private:
		float Q;
		float R;
	};

	class ObjectPlane : public Volume<float>
	{

	public:

		ObjectPlane(Vector x0, Vector n) : X0(x0), N(n.unitvector()) {}
		~ObjectPlane() {};

		const float eval(const Vector& X)const
		{

			float final = -1.0 * N * (X - X0);
			return final;
		}

	private:
		Vector X0;
		Vector N;
	};

	class ObjectCone : public Volume<float>
	{

	public:

		ObjectCone(float h, float thi_max, Vector n) : H(h), thi(thi_max), N(n) {}
		~ObjectCone() {};

		const float eval(const Vector& X)const
		{
			float temp = X * N;
			if (temp < 0)
			{
				return temp;
			}
			else if (temp > H)
			{
				return H - temp;
			}
			else if (0 < temp&&temp < H)
			{
				return X * N - X.magnitude() * cos(thi);
			}
			//return 0;
		}

	private:
		float thi, H;
		Vector N;
	};

	class ObjectCylinder : public Volume<float>
	{

	public:

		ObjectCylinder(float R, Vector N) : r(R), n(N) {}
		~ObjectCylinder() {};

		const float eval(const Vector& x)const
		{

			float final = r - (x - (x * n) * n).magnitude();
			return final;
		}

	private:
		float r;
		Vector n;
	};

	class ObjectIcosahedron : public Volume<float>
	{

	public:

		ObjectIcosahedron(float T) : t(T) {}
		~ObjectIcosahedron() {};

		const float eval(const Vector& x)const
		{
			float temp = x.magnitude();
			if (temp <= 1.8 * PI)
			{
				float final = cos(x.X() + t * x.Y()) + cos(x.X() - t * x.Y())
							+ cos(x.Y() + t * x.Z()) + cos(x.Y() - t * x.Z())
							+ cos(x.Z() - t * x.X()) + cos(x.Z() + t * x.X()) -2;
							return final;
			}

			else if (temp > 1.8 * PI)
			{
				return -1.8*PI;
			}
		}

	private:
		float t;
	};

	class ObjectSteinerPatch : public Volume<float>
	{

	public:

		ObjectSteinerPatch()  {}
		~ObjectSteinerPatch() {};

		const float eval(const Vector& x)const
		{

			float final = -1.0 * (x.X() * x.X() * x.Y() * x.Y() + x.X() * x.X() * x.Z() * x.Z()
						  +x.Y() * x.Y() * x.Z() * x.Z() - x.X() * x.Y() * x.Z());
			return final;
		}

	private:

	};

//==================================================================================================


//==============================================================================

}


#endif

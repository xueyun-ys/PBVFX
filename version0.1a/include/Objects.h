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
	class ObjectSFTranslate : public Volume<float>
	{

	public:

		ObjectSFTranslate(Volume<float> *f, Vector dx) : elem(f), dx(dx) {}
		~ObjectSFTranslate() {};

		const float eval(const Vector& x)const
		{
			//std::cout<<"Vector1:    "<<x.X()<<" "<<x.Y()<<" "<<x.Z()<<std::endl;
			//std::cout<<"Vector2:    "<<(x-dx).X()<<" "<<(x-dx).Y()<<" "<<(x-dx).Z()<<std::endl;
			return elem->eval(x - dx);
		}

	private:
		Volume<float> *elem;
		Vector dx;
	};

	class ObjectSFRotation : public Volume<float>
	{

	public:

		ObjectSFRotation(Volume<float>* f, Vector axis, float angle) : elem(f), ax(axis.unitvector()), a(angle) {}//normalization or unitvector
		~ObjectSFRotation() {};

		const float eval(const Vector& p)const
		{
			Vector vec_rotated = p * cos(a) + ax * (ax * p) * (1.0 - cos(a)) + (p ^ ax) * sin(a);
			float final = elem->eval(vec_rotated);
			return final;
		}

	private:
		Volume<float>* elem;
		Vector ax;
		float a;
	};

	class ObjectSFScale : public Volume<float>
	{

	public:

		ObjectSFScale(Volume<float>* f, Vector pivot, float scale) : elem(f), p(pivot), s(scale) {}
		~ObjectSFScale() {};

		const float eval(const Vector& x) const
		{
			return elem->eval((x - p)/s+p);
		}

	private:
		Volume<float>* elem;
		Vector p;
		float s;
	};
	class ObjectSFScale2 : public Volume<float>
	{

	public:

		ObjectSFScale2(Volume<float>* f, float scale, int type) : elem(f), s(scale), type(type) {}
		~ObjectSFScale2() {};

		const float eval(const Vector& x) const
		{
			if(type == 0)
			{
				Vector temp(x.X()/s, x.Y(), x.Z());
				return elem->eval(temp);
			}

			else if(type == 1)
			{
				Vector temp(x.X(), x.Y()/s, x.Z());
				return elem->eval(temp);
			}

			else if(type == 2)
			{
				Vector temp(x.X(), x.Y(), x.Z()/s);
				return elem->eval(temp);
			}

		}

	private:
		Volume<float>* elem;
		int type;
		float s;
	};
	//f(  f(   f(x)    )  )

	//=======================================================================================
	//operator redefinition......should build another file
	//color * scalar
	class ObjectColorDensityFieldMultiple : public Volume<Color>//return type is whatever in the <>
	{

	public:

		ObjectColorDensityFieldMultiple(Volume<float>* f, Volume<Color>* c) : elem(f), color(c) {}
		~ObjectColorDensityFieldMultiple() {};

		const Color eval(const Vector& x)const
		{
			//std::cout<<"hello11111111111111!"<<std::endl;

			 return Color(elem->eval(x) * color->eval(x).X(), elem->eval(x) * color->eval(x).Y(),
			 					elem->eval(x) * color->eval(x).Z(), elem->eval(x) * color->eval(x).W());
			//return Color(0.0, 0.0, 0.5,0.5);
		}

	private:
		Volume<float>* elem;
		Volume<Color>* color;
	};

	// color
	class ConstantColor : public Volume<Color>
	{
	public:
		ConstantColor(Color color) : constantC(color) {}
		~ConstantColor() {}

		const Color eval(const Vector& x) const { /*std::cout<<"hello2222222222222222222!"<<std::endl;*/ return constantC; }

	private:
		Color constantC;// = new Color(0,0,0.5,0);
	};

	//color+color
	class ObjectColorFieldAdd : public Volume<Color>//return type is whatever in the <>
	{

	public:

		ObjectColorFieldAdd(Volume<Color>* f, Volume<Color>* c) : elem(f), color(c) {}
		~ObjectColorFieldAdd() {};

		const Color eval(const Vector& x)const
		{
			return Color(elem->eval(x).X() + color->eval(x).X(), elem->eval(x).Y() + color->eval(x).Y(),
								elem->eval(x).Z() + color->eval(x).Z(), elem->eval(x).W() + color->eval(x).W());
		}

	private:
		Volume<Color>* elem;
		Volume<Color>* color;
	};

	//scalar + scalar
	class ObjectDensityFieldAdd : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectDensityFieldAdd(Volume<float>* f, Volume<float>* density) : elem(f), dens(density) {}
		~ObjectDensityFieldAdd() {};

		const float eval(const Vector& x)const
		{
			return (elem->eval(x) + dens->eval(x));
		}

	private:
		Volume<float>* elem;
		Volume<float>* dens;
	};
	//density_mask

//==============================================================================
	//Geo operator
	//addition
	class ObjectCombined : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectCombined(Volume<float>* a, Volume<float>* b) : elem1(a), elem2(b) {}
		~ObjectCombined() {};

		const float eval(const Vector& x)const
		{
			return max(elem1->eval(x), elem2->eval(x));
		}

	private:
		Volume<float>* elem1;
		Volume<float>* elem2;
	};

	//cutout
	class ObjectCutout : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectCutout(Volume<float>* a, Volume<float>* b) : elem1(a), elem2(b) {}
		~ObjectCutout() {};

		const float eval(const Vector& x)const
		{
			return min(elem1->eval(x), elem2->eval(x)*-1.0f);
		}

	private:
		Volume<float>* elem1;
		Volume<float>* elem2;
	};

	//intersection
	class ObjectIntersection : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectIntersection(Volume<float>* a, Volume<float>* b) : elem1(a), elem2(b) {}
		~ObjectIntersection() {};

		const float eval(const Vector& x)const
		{
			return min(elem1->eval(x), elem2->eval(x));
		}

	private:
		Volume<float>* elem1;
		Volume<float>* elem2;
	};

	//shell
	class ObjectShell : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectShell(Volume<float>* a, float c) : elem1(a), shell(c) {}
		~ObjectShell() {};

		const float eval(const Vector& x)const
		{
			return min(elem1->eval(x) + shell / 2.0f, shell / 2.0f - elem1->eval(x));
		}

	private:
		Volume<float>* elem1;
		float shell;
	};

	//Dilation
	class ObjectDilation : public Volume<float>//return type is whatever in the <>
	{

	public:

		ObjectDilation(Volume<float>* a, float c) : elem1(a), dil(c) {}
		~ObjectDilation() {};

		float eval(Vector& x)
		{
			return (elem1->eval(x) + dil);
		}

	private:
		Volume<float>* elem1;
		float dil;
	};

	//blend

//============================================================================
	//mask
	class ObjectMask : public Volume<float>//return type is whatever in the <>
	{
	public:

		ObjectMask(Volume<float>* a) : elem1(a) {}
		~ObjectMask() {};

		const float eval (const Vector& x) const
		{
			/*std::cout<<"hello333333333333!"<<std::endl;*/
			if (elem1->eval(x) > 0)
			{
				//std::cout<<"mask:    1"<<std::endl;
				return 1.0;
			}
			else if (elem1->eval(x) <= 0)
			{
				//std::cout<<"mask:    0"<<std::endl;
				return 0.0;
			}
		}

	private:
		Volume<float>* elem1;
	};

	//clamp
	class ObjectClamp : public Volume<float>//return type is whatever in the <>
	{
	public:

		ObjectClamp(Volume<float>* f, float a, float b) : elem1(f), a(a), b(b) {}
		~ObjectClamp() {};

		float eval(Vector& x)
		{
			if (elem1->eval(x) < a)
			{
				return a;
			}
			else if (elem1->eval(x) > b)
			{
				return b;
			}
			else if ( elem1->eval(x) >= a && elem1->eval(x) <= b)
			{
				return elem1->eval(x);
			}
		}

	private:
		Volume<float>* elem1;
		float a, b;//boundary
	};

	//========================================================================================


}


#endif

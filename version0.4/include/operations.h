#pragma once
#ifndef ____PBE_OPERATIONS_H____
#define ____PBE_OPERATIONS_H____

#include "Volume.h"
#include <cmath>


using namespace std;

namespace lux
{
//======================================transformation========================================================
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

	//==================================operator redefinition===========================================
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

	//color * density(scalar)
	//return elem->eval(x)*color->eval(x)
	class ObjectColorDensityFieldMultiple : public Volume<Color>//return type is whatever in the <>
	{

	public:

		ObjectColorDensityFieldMultiple(Volume<float>* f, Volume<Color>* c) : elem(f), color(c) {}
		~ObjectColorDensityFieldMultiple() {};

		const Color eval(const Vector& x)const
		{
			//std::cout<<"hello11111111111111!"<<std::endl;
			float temp = elem->eval(x);
			Color c(color->eval(x));
			 return Color(temp * c.X(), temp * c.Y(),
			 					temp * c.Z(), temp * c.W());
			//return Color(0.0, 0.0, 0.5,0.5);
		}

	private:
		Volume<float>* elem;
		Volume<Color>* color;
	};

	//color+color
	class ObjectColorFieldAdd : public Volume<Color>//return type is whatever in the <>
	{

	public:

		ObjectColorFieldAdd(Volume<Color>* f, Volume<Color>* c) : elem(f), color(c) {}
		~ObjectColorFieldAdd() {};

		const Color eval(const Vector& x)const
		{
			Color c(color->eval(x));
			Color e(elem->eval(x));
			return Color(e.X() + c.X(), e.Y() + c.Y(),
								e.Z() + c.Z(), e.W() + c.W());
		}

	private:
		Volume<Color>* elem;
		Volume<Color>* color;
	};

	//scalar + scalar/Density + Density
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


//==============================Geo operator=======================================
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
	class ObjectIntersection : public Volume<float>
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
	class ObjectShell : public Volume<float>
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
	class ObjectDilation : public Volume<float>
	{

	public:

		ObjectDilation(Volume<float>* a, float c) : elem1(a), dil(c) {}
		~ObjectDilation() {};

		const float eval(const Vector& x)const
		{
			return (elem1->eval(x) + dil);
		}

	private:
		Volume<float>* elem1;
		float dil;
	};

	//blend

//============================mask===========================================
	//density_mask
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
			else// if (elem1->eval(x) <= 0)
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

		const float eval(const Vector& x)const
		{//one eval
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

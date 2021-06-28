//*******************************************************************
//
//   CreateFields.h
//
//
//
//
//
//*******************************************************************

#ifndef __CREATEFIELDS_H__
#define __CREATEFIELDS_H__

#include <cmath>
#include <string>
#include <iostream>

#include "Objects.h"
#include "operations.h"
#include "Color.h"
#include "Image.h"


using namespace img;
using namespace std;
using namespace lux;

namespace lux
{

class Avator
{
  public:

   Avator()
   {
     //create objects
     //head
         //ConstantField *densityVolume = new ConstantField(1.0);
         //eyes
       	ObjectSphere *densityVolume2 = new ObjectSphere(2.5);
         ObjectEllipse *densityVolume = new ObjectEllipse(Vector(1.0, 0.0, 0.0), 5.0, 2.0 );
         ObjectIntersection *s4 = new ObjectIntersection(densityVolume2, densityVolume);
         ObjectSFScale *s5 = new ObjectSFScale(s4, Vector(0.0, 0.0, 0.0), 0.3);
         ObjectSFTranslate *s6 = new ObjectSFTranslate(s5, Vector(2.0, 0.0, 0.0));
         ObjectCombined *s7 = new ObjectCombined(s5, s6);
         ObjectSFTranslate *t_eyes = new ObjectSFTranslate(s7, Vector(0.0, 3.5, 3.0));
         //head
         ObjectSphere *d3 = new ObjectSphere(2.5);
         ObjectSFTranslate *s1 = new ObjectSFTranslate(d3, Vector(1.0, 3.3, 1.0));
         //head+eyes
         ObjectMask *m2 = new ObjectMask(t_eyes);
         ObjectMask *m3 = new ObjectMask(s1);
         ObjectDensityFieldAdd *m1 = new ObjectDensityFieldAdd(m2, m3);

         Color c1(0.5d,0.3d,0.0d,0.0d);
         Color c2(0.2d,0.2d,0.4d,0.0d);
       	 ConstantColor *colorVolume = new ConstantColor(c1);
         ConstantColor *colorVolume2 = new ConstantColor(c2);
         //colorVolume->set(0.0, 0.0, 0.5, 0.0);

         ObjectColorDensityFieldMultiple *cf1 = new ObjectColorDensityFieldMultiple(m2, colorVolume);
         ObjectColorDensityFieldMultiple *cf2 = new ObjectColorDensityFieldMultiple(m3, colorVolume2);
         ObjectColorFieldAdd *s2 = new ObjectColorFieldAdd(cf1, cf2);
     //head end
     //hair_ponytail?
         ObjectCone *d4 = new ObjectCone(6.0, 0.2, Vector(0.0, -1.0, 0.0));
         ObjectSFRotation *s8 = new ObjectSFRotation(d4, Vector(0.0, 0.0, 1.0), 3.1415*2.0/3.0);
         ObjectSFTranslate *t1 = new ObjectSFTranslate(s8, Vector(4.0, 1.7, 0.0));

         ObjectBox *d5 = new ObjectBox(1.0, 5.0);
         ObjectSFRotation *t2 = new ObjectSFRotation(d5, Vector(0.0, 0.0, 1.0), 3.1415/6.0);
         ObjectSFTranslate *t3 = new ObjectSFTranslate(t2, Vector(5.0, 2.5, 0.0));
         ObjectCombined *s9 = new ObjectCombined(t3, t1);

         ObjectCone *s10 = new ObjectCone(3.0, 0.4, Vector(0.0, -1.0, 0.0));
         ObjectSFRotation *s11 = new ObjectSFRotation(s10, Vector(0.0, 0.0, 1.0), -3.1415/3.0);
         ObjectSFScale *s12 = new ObjectSFScale(s11, Vector(0.0, 0.0, 0.0), 0.7);
         ObjectSFTranslate *t4 = new ObjectSFTranslate(s12, Vector(5.0, 2.5, 0.0));

         ObjectCombined *s13 = new ObjectCombined(s9, t4);
         ObjectSFScale *s14 = new ObjectSFScale(s13, Vector(0.0, 0.0, 0.0), 0.4);

         ObjectTorus *s15 = new ObjectTorus(Vector(0.0, 0.0, 1.0), 5.0, 1.0);
         ObjectPlane *s16 = new ObjectPlane(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0));
         ObjectCutout *s17 = new ObjectCutout(s15, s16);
         ObjectPlane *s18 = new ObjectPlane(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0));
         ObjectSFTranslate *s19 = new ObjectSFTranslate(s18, Vector(0.0, -2.0, 0.0));
         ObjectCutout *s20 = new ObjectCutout(s17, s19);
         ObjectSFScale *s21 = new ObjectSFScale(s20, Vector(0.0, 0.0, 0.0), 0.27);
         ObjectSFRotation *t5 = new ObjectSFRotation(s21, Vector(0.0, 0.0, 1.0), 3.1415/6.0);

         ObjectSFTranslate *s22 = new ObjectSFTranslate(t5, Vector(4.2, 0.6, -0.2));
         ObjectSFTranslate *s23 = new ObjectSFTranslate(t5, Vector(4.2, 0.6, 0.2));
         ObjectCombined *s24 = new ObjectCombined(s22, s23);

         ObjectCombined *s26 = new ObjectCombined(s14, s24);
         ObjectSFTranslate *s27 = new ObjectSFTranslate(s26, Vector(0.0, 4.0, 0.0));
         ObjectSFRotation *s28 = new ObjectSFRotation(s27, Vector(0.0, 1.0, 0.0), 3.1415/2.0);
         ObjectMask *m4 = new ObjectMask(s28);

         Color c3(0.2d,0.6d,0.3d,0.0d);
         ConstantColor *colorVolume3 = new ConstantColor(c3);
         ObjectColorDensityFieldMultiple *cf3 = new ObjectColorDensityFieldMultiple(m4, colorVolume3);

         ObjectColorFieldAdd *color1 = new ObjectColorFieldAdd(cf3, s2);

         ObjectDensityFieldAdd *m5 = new ObjectDensityFieldAdd(m1, m4);


         //body
         ObjectBox *s29 = new ObjectBox(1.0, 5.0);
         ObjectSFScale2 *s30 = new ObjectSFScale2(s29, 0.7, 2);
         ObjectSFScale *s31 = new ObjectSFScale(s30, Vector(0.0, 0.0, 0.0), 3.0);
         ObjectSFTranslate *s32 = new ObjectSFTranslate(s31, Vector(1.0, -1.3, 1.0));
         ObjectCutout *s33 = new ObjectCutout(s32, s1);

         ObjectMask *m6 = new ObjectMask(s33);
         Color c4(0.4d,0.6d,0.7d,0.0d);
         ConstantColor *colorVolume4 = new ConstantColor(c4);
         ObjectColorDensityFieldMultiple *cf4 = new ObjectColorDensityFieldMultiple(m6, colorVolume4);

         ObjectColorFieldAdd *color2 = new ObjectColorFieldAdd(cf4, color1);
         ObjectDensityFieldAdd *m7 = new ObjectDensityFieldAdd(m5, m6);

     //feet
         ObjectIcosahedron *s34 = new ObjectIcosahedron(1.0);
         ObjectSFScale *s35 = new ObjectSFScale(s34, Vector(0.0, 0.0, 0.0), 0.2);
         ObjectSFTranslate *s36 = new ObjectSFTranslate(s35, Vector(0.0, -5.0, 0.0));
         ObjectSFTranslate *s37 = new ObjectSFTranslate(s35, Vector(2.8, -5.0, 0.0));
         ObjectCombined *s38 = new ObjectCombined(s36, s37);
         //ObjectCutout *s39 = new ObjectCutout(s38, s32);

         ObjectMask *m8 = new ObjectMask(s38);
         Color c5(0.4d,0.8d,0.7d,0.0d);
         ConstantColor *colorVolume5 = new ConstantColor(c5);
         ObjectColorDensityFieldMultiple *cf5 = new ObjectColorDensityFieldMultiple(m8, colorVolume5);

         ObjectColorFieldAdd *color3 = new ObjectColorFieldAdd(cf5, color2);
         ObjectDensityFieldAdd *m9 = new ObjectDensityFieldAdd(m7, m8);

     //arms
         ObjectEllipse *s40 = new ObjectEllipse(Vector(1.0, 0.0, 0.0), 5.0, 2.0 );
         ObjectSFScale *s41 = new ObjectSFScale(s40, Vector(0.0, 0.0, 0.0), 0.5);
         ObjectSFRotation *ts44 = new ObjectSFRotation(s41, Vector(0.0, 0.0, 1.0), -3.1415/3.0);
         ObjectSFRotation *s43 = new ObjectSFRotation(s41, Vector(0.0, 0.0, 1.0), 3.1415/3.0);
         ObjectSFTranslate *s44 = new ObjectSFTranslate(s43, Vector(5.3, 0.0, 0.0));
         ObjectSFTranslate *s42 = new ObjectSFTranslate(ts44, Vector(-3.5, 0.0, 0.0));

         ObjectCombined *s45 = new ObjectCombined(s42, s44);
         //ObjectCutout *s39 = new ObjectCutout(s38, s32);

         ObjectMask *m10 = new ObjectMask(s45);
         Color c6(0.6d,0.5d,0.7d,0.0d);
         ConstantColor *colorVolume6 = new ConstantColor(c6);
         ObjectColorDensityFieldMultiple *cf6 = new ObjectColorDensityFieldMultiple(m10, colorVolume6);

         ObjectColorFieldAdd *color4 = new ObjectColorFieldAdd(cf6, color3);
         ObjectDensityFieldAdd *m11 = new ObjectDensityFieldAdd(m9, m10);

     //weapon
         ObjectCylinder *s46 = new ObjectCylinder(0.5, Vector(0.0, 1.0, 0.0));
         ObjectSFTranslate *s47 = new ObjectSFTranslate(s46, Vector(-4.5, 0.0, 0.0));
         //ObjectSphere *s48 = new ObjectSphere(20.0);
         //ObjectIntersection *s49 = new ObjectIntersection(s47, s48);
         ObjectPlane *s48 = new ObjectPlane(Vector(0.0, 4.5, 0.0), Vector(0.0, -1.0, 0.0));
         ObjectPlane *s49 = new ObjectPlane(Vector(0.0, -6.0, 0.0), Vector(0.0, 1.0, 0.0));
         ObjectCutout *s50 = new ObjectCutout(s47, s48);
         ObjectCutout *s51 = new ObjectCutout(s50, s49);

         ObjectMask *m12 = new ObjectMask(s51);
         Color c7(0.6d,0.7d,0.2d,0.0d);
         ConstantColor *colorVolume7 = new ConstantColor(c7);
         ObjectColorDensityFieldMultiple *cf7 = new ObjectColorDensityFieldMultiple(m12, colorVolume7);

         ObjectColorFieldAdd *color5 = new ObjectColorFieldAdd(cf7, color4);
         ObjectDensityFieldAdd *m13 = new ObjectDensityFieldAdd(m11, m12);


         ObjectSteinerPatch *s52 = new ObjectSteinerPatch();
         ObjectSFScale *ss1 = new ObjectSFScale(s52, Vector(0.0, 0.0, 0.0), 3.0);
         ObjectSFTranslate *s53 = new ObjectSFTranslate(ss1, Vector(-4.5, 6.0, 0.0));


         ObjectMask *m14 = new ObjectMask(s53);
         Color c8(0.5d,0.2d,0.6d,0.0d);
         ConstantColor *colorVolume8 = new ConstantColor(c8);
         ObjectColorDensityFieldMultiple *cf8 = new ObjectColorDensityFieldMultiple(m14, colorVolume8);

         color6 = new ObjectColorFieldAdd(cf8, color5);
         m15 = new ObjectDensityFieldAdd(m13, m14);

         //ObjectShell *s4 = new ObjectShell(densityVolume2, 0.2);
   }

   ~Avator(){}

   ObjectDensityFieldAdd* getfield()
   {
      return m15;
   }
   ObjectColorFieldAdd* getcolor()
   {
      return color6;
   }

  private:
    ObjectDensityFieldAdd *m15;
    ObjectColorFieldAdd *color6;
};

}



#endif

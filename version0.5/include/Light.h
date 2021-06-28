//*******************************************************************
//
//   Light.h
//
// Light class in the namespace lux
//
//
//
//*******************************************************************

#ifndef __LUX_LIGHT_H__
#define __LUX_LIGHT_H__

#include "grids.h"

namespace lux
{

class Light
{
  public:

   Light(Color C, Vector P) : c(C), p(P)//, gtf(GTF), cf(color_field)
   {

   }

   ~Light(){}

   Vector& pos() { return p; }
   Color& color() {return c;}
   //GridsToFields* get_dense() {return gtf;}
   //ObjectColorFieldAdd* get_cf() {return cf;}

   private:
   Color c;
   Vector p;
   //GridsToFields *gtf;
   //ObjectColorFieldAdd * cf;
};

// void CreateLights()
// {
//   Grids *gr_light = new Grids(Vector(-6.0, -6.0, 6.0), 61,61,61,0.2);
//   DSM_Grids *dsm_grids = new  DSM_Grids(*gr_light, l, gtf);
// }

}



#endif

//*******************************************************************
//
//   grids.h
//
// grids for speed up, at the cost of quality
//
//
//
//*******************************************************************
#pragma once
#ifndef ____GRIDS_H____
#define ____GRIDS_H____

#include <stdlib.h>
#include "Light.h"
//#include "Volume.h"
using namespace std;

namespace lux
{

    //float grids
    class Grids
    {
      public:

       Grids(Vector llc, int x, int y, int z, float step) : llc(llc)
       {
         xyz[0] = x;
         xyz[1] = y;
         xyz[2] = z;
         delta = step;
         size[0] = (float)(x-1) * step;
         size[1] = (float)(y-1) * step;
         size[2] = (float)(z-1) * step;
         data.resize(x*y*z);
       }
       Grids()
       {
         xyz[0] = xyz[1] = xyz[2] = 20;
       }
       ~Grids(){}

       //!  Set all four components
       void set( const float vx, const float vy, const float vz)
       {
          xyz[0] = vx;
          xyz[1] = vy;
          xyz[2] = vz;
       }

       Vector pos(int i, int j, int k)
       {
         return Vector(llc.X() + i * delta, llc.Y() + j * delta, llc.Z() - k * delta);
       }

       float trilinearInterpolation(Vector P)
       {
         int i, j, k;
         //int i = int((P.X()-llc.X())/delta);
         float ii = P.X()-llc.X();
         float jj = P.Y()-llc.Y();
         float kk = llc.Z()-P.Z();
         if(ii<0||jj<0||kk<0)
         {
           return 0.0;
         }
         else
         {
           i = (int)(ii/delta);
           j = (int)(jj/delta);
           k = (int)(kk/delta);
         }
         if(i<xyz[0] && j<xyz[1] && k<xyz[2])
         {
           float wi = (float)(P.X()-llc.X())/delta - (float)i;
           float wj = (float)(P.Y()-llc.Y())/delta - (float)j;
           float wk = (float)(llc.Z()-P.Z())/delta - (float)k;

           //std::cout<<"wi: "<<wi<<"     wj: "<<wj<<"      wk: "<<wk<<std::endl;
           //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

           float result = gridValue(i, j, k)*(1-wi)*(1-wj)*(1-wk)+
           gridValue(i+1, j, k)*wi*(1-wj)*(1-wk)+
           gridValue(i, j+1, k)*(1-wi)*wj*(1-wk)+
           gridValue(i, j, k+1)*(1-wi)*(1-wj)*wk+
           gridValue(i+1, j+1, k)*wi*wj*(1-wk)+
           gridValue(i+1, j, k+1)*wi*(1-wj)*wk+
           gridValue(i, j+1, k+1)*(1-wi)*wj*wk+
           gridValue(i+1, j+1, k+1)*wi*wj*wk;

           return result;
         }
         else
         {
           return -0.0;//can't set to zero, zero = CSG edge
         }
       }

       float gridValue(int i, int j, int k)
       {
         return data[i+j*xyz[0]+k*xyz[1]*xyz[0]];
       }

       float gridindex(int i, int j, int k)
       {
         return i+j*xyz[0]+k*xyz[1]*xyz[0];
       }

       const int X() const { return xyz[0]; }
       const int Y() const { return xyz[1]; }
       const int Z() const { return xyz[2]; }
       float get_delta() const {return delta;}
       Vector get_llc() const {return llc;}
       std::vector<float>& get_data() {return data;}//const

      private:

      Vector llc;
      int xyz[3];
      float delta;
      float size[3];
      //float* data*;//???
      std::vector<float> data;
    };


    class FieldsToGrids
    {
       public:

       FieldsToGrids(Grids G, Volume<float> *F) : g(G), f(F)
       {
          float s = g.get_delta();
          int test = 0;
        	for(int z=0; z< g.Z(); z++)
        	{
          	  for(int y = 0; y < g.Y(); y++)
          	  {
          	    for(int x = 0; x< g.X(); x++)
          	    {
          	      float n = g.gridindex(x, y, z);
          	      Vector p = g.pos(x, y, z);
          	      g.get_data()[n] = f->eval(p);
          	    }
          	  }
        	 }
       }
       //FieldsToGrids();
       // {
       //   ;
       // }
       ~FieldsToGrids(){}

       Grids get_grid()
       {
         return g;
       }

       private:

       Grids g;
       Volume<float> *f;
    };

    class GridsToFields
    {
      public:

       GridsToFields(Grids G) : g(G)
       {

       }
       ~GridsToFields(){}

       float eval(Vector P) //const//wrong code
       {
         float res=g.trilinearInterpolation(P);
         return res;
       }

       private:

       Grids g;
       //DSM_Grids dsm;
    };

    class DSM_Grids
    {
       public:

       DSM_Grids(Grids G, Light L, GridsToFields *gtf) : g(G), l(L), gtf(gtf)
       {

       }

       ~DSM_Grids(){}

       float dsm_ijk(Vector p)
       {
         GridsToFields *dens = gtf;
         if(dens->eval(p)>0)
         {
           //light ray march pos
           Vector X = p;
           //light position
           Vector lp = l.pos();
           //light direction
           Vector temp = lp - X;
           float distance = temp.magnitude();
           Vector nl = temp.unitvector();
           //float Tl = 1.0;
           //float Kl = 0.9;
           double del_sl = 0.05;//0.01
           float dsm = 0.0;
           //float alpha;
           double s = 0.0;
           //Lp = (0, 0, 0, 0)
           //Color Lp;

           //ObjectColorFieldAdd *color_field = l.get_cf();
           while (s < distance)
           {
             X = X + del_sl * nl;//dei = 0 is missing;
             s += del_sl;
             dsm += del_sl * dens->eval(X);
             //del_t = exp(-1.0 * K * del_s * gtf->eval(X));
             //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
             //Lp = Lp + color_field->eval(X) * T * (1.0 - del_t) / K;
             //T *= del_t;
             //Lp.set(Lp.X(), Lp.Y(), Lp.Z(), alpha);
           }
           //alpha = 1.0 - T;
           return dsm;
         }
         return 0.0;
       }
    //float a*;//wrong
       void bakeDSM()
       {
         float s = g.get_delta();
         int test = 0;
         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               float n = g.gridindex(x, y, z);
               Vector p = g.pos(x, y, z);
               float temp = dsm_ijk(p);
               g.get_data()[n] = temp;
               //g.get_data()[n] = dsm_ijk(p);

               //test++;
               //std::cout<<"step: "<<test<<"  "<<temp<<endl;
             }
           }
          }
        }

        Grids get_grid()
        {
          return g;
        }

       private:

       GridsToFields *gtf;
       Grids g;
       //Vector p;
       Light l;
    };

    class GridsToLightFields
    {
      public:

       //GridsToLightFields(DSM_Grids DSM) : dsm(DSM)
       GridsToLightFields(Grids g) : g(g)
       {

       }
       ~GridsToLightFields(){}

       float light_eval(Vector P) //const//wrong code
       {
         //Grids g = dsm.get_grid();//segment fault(core dumped)
         float res = g.trilinearInterpolation(P);
         //kappa
         res = exp(-0.9*res);
         //std::cout<<"step: "<<"  "<<res<<endl;
         return res;
       }

       private:

       Grids g;
       //DSM_Grids dsm;//must have parameters at the start of this construction function
    };

}
#endif

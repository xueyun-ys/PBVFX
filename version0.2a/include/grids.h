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
         //std::cout<<"posx:"<<llc.X() + i * delta<<" posy "<<llc.Y() + j * delta<<" posz "<<llc.Z() - k * delta<<std::endl;
         return Vector(llc.X() + i * delta, llc.Y() + j * delta, llc.Z() - k * delta);
       }

       float trilinearInterpolation(Vector P)
       {
         //std::cout<<" Px: "<<P.X()<<"  Py  "<<P.Y()<<"  Pz  "<<P.Z()<<std::endl;
         int i, j, k;
         //int i = int((P.X()-llc.X())/delta);
         float ii = P.X()-llc.X();
         float jj = P.Y()-llc.Y();
         float kk = llc.Z()-P.Z();
         //std::cout<<"22222222      jj:"<<jj<<std::endl;
         if(ii<0||jj<0||kk<0)
         {
           //std::cout<<"11111111"<<std::endl;
           return 0.0;
         }
         else
         {
           //std::cout<<"22222222      jj:"<<j<<std::endl;
           i = (int)(ii/delta);
           j = (int)(jj/delta);
           k = (int)(kk/delta);
           //std::cout<<"22222222      j:"<<j<<std::endl;
         }
//std::cout<<"test1111"<<std::endl;
//std::cout<<"i: "<<i<<"      j: "<<j<<"      k: "<<k<<std::endl;
         if(i<xyz[0] && j<xyz[1] && k<xyz[2])
         {
           float wi = (float)(P.X()-llc.X())/delta - (float)i;
           float wj = (float)(P.Y()-llc.Y())/delta - (float)j;
           float wk = (float)(llc.Z()-P.Z())/delta - (float)k;
           //std::cout<<" Py: "<<P.Y()<<"  llcy  "<<llc.Y()<<std::endl;
           //std::cout<<"wj: "<<(float)(P.Y()-llc.Y())/delta<<"     j: "<<j<<std::endl;
           //std::cout<<"wi: "<<wi<<"     wj: "<<wj<<"      wk: "<<wk<<std::endl;
           //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
           // if(wi==0&&wj==0&&wk==0)
           // {
           //   return gridValue(i, j, k);
           // }
           //system("pause");//windows
           float result = gridValue(i, j, k)*(1-wi)*(1-wj)*(1-wk)+
           gridValue(i+1, j, k)*wi*(1-wj)*(1-wk)+
           gridValue(i, j+1, k)*(1-wi)*wj*(1-wk)+
           gridValue(i, j, k+1)*(1-wi)*(1-wj)*wk+
           gridValue(i+1, j+1, k)*wi*wj*(1-wk)+
           gridValue(i+1, j, k+1)*wi*(1-wj)*wk+
           gridValue(i, j+1, k+1)*(1-wi)*wj*wk+
           gridValue(i+1, j+1, k+1)*wi*wj*wk;
           //std::cout<<"gridvalue!!!"<<gridValue(i+1, j, k)<<std::endl;
           //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
           //std::cout<<"!!!!!!!!!"<<result<<std::endl;
           //if(result>0){std::cout<<"nice!!!"<<std::endl;}
           return result;
         }
         else
         {
           //std::cout<<"======="<<std::endl;

           return -0.0;//can't set to zero, zero = CSG edge
         }
       }

       float gridValue(int i, int j, int k)
       {
         //std::cout<<"x: "<<xyz[0]<<" y: "<<xyz[1]<<std::endl;
         //std::cout<<"i: "<<i<<" j: "<<j<<"k:"<<k<<std::endl;
         //std::cout<<i+j*xyz[0]+k*xyz[1]*xyz[0]<<std::endl;
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
      //float data*;
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
                  //test++;
                  //std::cout<<"x: "<<x<<"y: "<<y<<"z: "<<z<<std::endl;
          	      float n = g.gridindex(x, y, z);
                  //std::cout<<"x: "<<x<<"y:  "<<y<<"z:  "<<z<<std::endl;
          	      Vector p = g.pos(x, y, z);//(g.get_llc().X()+ x*s, g.get_llc().Y()+y*s, g.get_llc().Z()-z*s);
                  //std::cout<<"x: "<<p.X()<<"      y: "<<p.Y()<<"      z: "<<p.Z()<<std::endl;
                  //Vector temp(1.3, 5.3, 1.8);
                  //std::cout<<"temp!!!!!: "<< f->eval(temp)<<std::endl;
                  //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
          	      g.get_data()[n] = f->eval(p);
                  //if(g.get_data()[n]!=0) {std::cout<<"x: "<<p.X()<<"      y: "<<p.Y()<<"      z: "<<p.Z()<<std::endl;
                  //if(f->eval(p)!=0) {std::cout<<"x: "<<p.X()<<"      y: "<<p.Y()<<"      z: "<<p.Z()<<"  "<<f->eval(p)<<"  "<<g.get_data()[n]<<std::endl;

                  //std::cout<< test<<":           "<<"n:"<<n<<"  "<<f->eval(p)<<std::endl;
                  //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');}
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
       GridsToFields()
       {
       }
       ~GridsToFields(){}

       const float eval(Vector P) //const//wrong code
       {
         //std::cout<<"test2222"<<std::endl;
         float res=g.trilinearInterpolation(P);
         return res;
       }

       private:

       Grids g;
    };
}
#endif

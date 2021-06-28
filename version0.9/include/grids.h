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
#include <algorithm>
#include "createOBJ.h"
#include "operations.h"
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

//==============================================================================
       float trilinearInterpolation(Vector P)
       {
//cout<<"hello11111111"<<endl;

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
//cout<<"hello2222222"<<endl;
         if(i<xyz[0]-1 && j<xyz[1]-1 && k<xyz[2]-1)
         {
           float wi = (float)(P.X()-llc.X())/delta - (float)i;
           float wj = (float)(P.Y()-llc.Y())/delta - (float)j;
           float wk = (float)(llc.Z()-P.Z())/delta - (float)k;
//cout<<"hello3333333"<<endl;
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
//cout<<"test0000000"<<endl;
           return result;
         }
         else
         {
           return -0.0;//can't set to zero, zero = CSG edge
         }
       }

       void BakeWisp(Vector P,  float rho)
       {
         int i, j, k;
         float ii = P.X()-llc.X();
         float jj = P.Y()-llc.Y();
         float kk = llc.Z()-P.Z();
         if(ii<0||jj<0||kk<0)
         {
           return;
         }
         else
         {
           i = (int)(ii/delta);
           j = (int)(jj/delta);
           k = (int)(kk/delta);
           // cout<<"P: "<<P.X()<<" "<<P.Y()<<" "<<P.Z()<<endl;
           // cout<<"i: "<<ii<<"j: "<<jj<<"k: "<<kk<<endl;
           // cout<<"ijk: "<<i<<" "<<j<<" "<<k<<" "<<endl;
         }
         if(i<xyz[0]-1 && j<xyz[1]-1 && k<xyz[2]-1)
         {
           float wi = (float)(ii)/delta - (float)i;
           float wj = (float)(jj)/delta - (float)j;
           float wk = (float)(kk)/delta - (float)k;
           //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

           //Bake into grids

           data[i+j*xyz[0]+k*xyz[1]*xyz[0]] += rho*(1.0-wi)*(1-wj)*(1-wk);
           data[i+1+j*xyz[0]+k*xyz[1]*xyz[0]] += rho*wi*(1-wj)*(1-wk);
           data[i+(j+1)*xyz[0]+k*xyz[1]*xyz[0]] += rho*(1.0-wi)*wj*(1-wk);
           data[i+j*xyz[0]+(k+1)*xyz[1]*xyz[0]] += rho*(1.0-wi)*(1-wj)*wk;
           data[(i+1)+(j+1)*xyz[0]+k*xyz[1]*xyz[0]] += rho*wi*wj*(1-wk);
           data[(i+1)+j*xyz[0]+(k+1)*xyz[1]*xyz[0]] += rho*wi*(1-wj)*wk;
           data[i+(j+1)*xyz[0]+(k+1)*xyz[1]*xyz[0]] += rho*(1.0-wi)*wj*wk;
           data[(i+1)+(j+1)*xyz[0]+(k+1)*xyz[1]*xyz[0]] += rho*wi*wj*wk;

           //cout<<"bake: "<< rho*(1.0-wi)*(1-wj)*(1-wk);
         }
       }
//==============================================================================
       float gridValue(int i, int j, int k)
       {
         //cout<<"hello3333333"<<"i "<<i<<"j "<<j<<"k "<<k<<endl;
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
       std::vector<float>& get_data() {return data;}//const &!!!
       //std::vector<double>& get_data() {return data;}//const &!!!

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

    class GridsToFields //: public Volume<float>
    {
      public:

       GridsToFields(Grids G) : g(G)
       {

       }
       ~GridsToFields(){}

       //const float eval(Vector P) const//wrong code
       float eval(Vector P)
       {
         ///Vector PP = P;
         float res=g.trilinearInterpolation(P);
         return res;
       }

       private:

       Grids g;
       //DSM_Grids dsm;
    };

//=================================lights==========================================
    class DSM_Grids
    {
       public:

       DSM_Grids(Grids G, Light L, GridsToFields *gtf) : g(G), l(L), gtf(gtf)
       //DSM_Grids(Grids G, Light L, Volume<float> *gtf) : g(G), l(L), gtf(gtf)
       {

       }
       // DSM_Grids(Grids G, Light L, ObjectDensityFieldAdd *dfa) : g(G), l(L), dfa(dfa)
       // {
       //
       // }

       ~DSM_Grids(){}

       float dsm_ijk(Vector p)
       {
         //GridsToFields *dens = gtf;
         if(gtf->eval(p)>0)//dens
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
             dsm += del_sl * gtf->eval(X);//dens
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
         //int test = 0;
         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               float n = g.gridindex(x, y, z);
               Vector p = g.pos(x, y, z);
               g.get_data()[n] = dsm_ijk(p);
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
       //Volume<float>* gtf;
       Grids g;
       //Vector p;
       Light l;
       //ObjectDensityFieldAdd *dfa;
    };
    class DSM_Grids2
    {
       public:
       DSM_Grids2(Grids G, Light L, Volume<float> *dfa) : g(G), l(L), dfa(dfa)
       {

       }

       ~DSM_Grids2(){}

       float dsm_ijk(Vector p)
       {
         //GridsToFields *dens = gtf;
         if(dfa->eval(p)>0)//dens
         {
           //light ray march pos
           Vector X = p;
           //light position
           Vector lp = l.pos();
           //light direction
           Vector temp = lp - X;
           float distance = temp.magnitude();
           Vector nl = temp.unitvector();

           double del_sl = 0.05;//0.01
           float dsm = 0.0;
           double s = 0.0;

           //ObjectColorFieldAdd *color_field = l.get_cf();
           while (s < distance)
           {
             X = X + del_sl * nl;//dei = 0 is missing;
             s += del_sl;
             dsm += del_sl * dfa->eval(X);//dens
           }
           return dsm;
         }
         return 0.0;
       }

       void bakeDSM()
       {
         float s = g.get_delta();
         //int test = 0;
         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               float n = g.gridindex(x, y, z);
               Vector p = g.pos(x, y, z);
               g.get_data()[n] = dsm_ijk(p);
             }
           }
         }
        }

        Grids get_grid()
        {
          return g;
        }

       private:

       Grids g;
       Light l;
       Volume<float> *dfa;
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
         //Grids g = dsm.get_grid();//segment fault(core dumped)???????
         float res = g.trilinearInterpolation(P);
         //kappa
         res = exp(Kappa*res);
         return res;
       }

       private:

       Grids g;
       float Kappa = -0.9;
       //DSM_Grids dsm;//must have parameters at the start of this construction function
    };

//==================================lights=========================================
//===================================levelset=========================================
    class LS_Grids
    {
       public:

       LS_Grids(OBJ_model obj, Grids G) : obj(obj), g(G)
       {
         //initialize grid
         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               float n = g.gridindex(x, y, z);
               g.get_data()[n] = 1000000.0;
             }
           }
         }

         vector<vector<int>> f = obj.getface();
         vector<vector<double>> vp = obj.getpoints();
         int size = f.size();
         double maxx, maxy, maxz;
         double minx, miny, minz;
         float delta = g.get_delta();
         //calculate shortest distance for each grid within bandwise
         int test = 0;
         //cout<<"size_f:  "<<size<<endl;
         //cout<<"size_vp:  "<<vp.size()<<endl;
         //int ptest = f[4961][0];

         //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
         //int point1, point2, point3;
         //double x1, y1, z1;
         double x2, y2, z2;
         //double x3, y3, z3;
         for(int i=0; i<size; i++)
         {
//           cout<<"test for levelset:  "<<test<<endl;
           test++;
           int point1, point2, point3;
           point1 = f[i][0]-1;
           point2 = f[i][1]-1;
           point3 = f[i][2]-1;
//cout<<"test444444444444:  "<<test<<endl;
           double x1, y1, z1;
           x1 = vp[point1][0];
           y1 = vp[point1][1];
           z1 = vp[point1][2];
//cout<<"test444444444444:  "<<vp[point2][0]<<" "<<vp[point2][1]<<" "<<vp[point2][2]<<endl;
//cout<<"point: "<<point2<<endl;
//cout<<"point: "<<point3<<endl;
           double x2, y2, z2;
           x2 = vp[point2][0];
           y2 = vp[point2][1];
           z2 = vp[point2][2];
//cout<<"test555555555555:  "<<endl;//<<vp[point3][0]<<" "<<vp[point3][1]<<" "<<vp[point3][2]<<endl;
           double x3, y3, z3;
           x3 = vp[point3][0];
           y3 = vp[point3][1];
           z3 = vp[point3][2];
//cout<<"test11111111111111:  "<<test<<endl;

           maxx = max(x1, x2);
           maxx = max(maxx, x3);
           maxy = max(y1, y2);
           maxy = max(maxy, y3);
           maxz = max(z1, z2);
           maxz = max(z3, maxz);
           minx = min(x1, x2);
           minx = min(minx, x3);
           miny = min(y1, y2);
           miny = min(miny, y3);
           minz = min(z1, z2);
           minz = min(minz, z3);
           // cout<<"test2222222222222:  zmin"<<minz<<"  zmax "<<maxz<<endl;
           // cout<<"test2222222222222:  ymin"<<miny<<"  ymax "<<maxy<<endl;
           // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

           // if(x1>x2)
           // {
           //   if(x1>x3)
           //   {
           //     maxx = x1;
           //   }
           //   else
           //   {
           //     maxx = x
           //   }
           // }
           //narrow band region
           int x_lower, x_upper, y_lower, y_upper, z_lower, z_upper;

           Vector t = g.get_llc();
           double temp = t.X();
           for(int j = 0; j< g.X(); j++)
           {
             //float n = g.gridindex(x, y, z);
             //Vector p = g.pos(x, y, z);
             if(minx>=temp + j*delta)
             {

             }
             else
             {
               x_lower = j-1;
               break;
             }
           }
           for(int j = 0; j< g.X(); j++)
           {
             //double temp = t.X();
             if(maxx<=temp + j*delta)
             {
               x_upper = j;
               break;
             }
             else
             {

             }
           }

           temp = t.Y();
           for(int j = 0; j< g.Y(); j++)
           {

             //double temp = t.Y();
             if(miny>=temp + j*delta)
             {

             }
             else
             {
               y_lower = j-1;
               break;
             }
           }
           for(int j = 0; j< g.Y(); j++)
           {
             //double temp = t.Y();
             if(maxy<=temp + j*delta)
             {
               y_upper = j;
               break;
             }
             else
             {

             }
           }

           temp = t.Z();
           for(int j = 0; j< g.Z(); j++)
           {

             //double temp = t.Z();
             if(minz < temp - j*delta)
             {

             }
             else
             {
               z_upper = j;
               break;
             }
           }
           for(int j = 0; j< g.Z(); j++)
           {
             //double temp = t.Z();
             if(maxz > temp - j*delta)
             {
               z_lower = j - 1;
               break;
             }
             else
             {

             }
           }
// cout<<"test2222222222222:  xlower"<<x_lower<<"  xupper "<<x_upper<<endl;
// cout<<"test2222222222222:  ylower"<<y_lower<<"  yupper "<<y_upper<<endl;
// cout<<"test2222222222222:  zlower"<<z_lower<<"  zupper "<<z_upper<<endl;
// do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

           //calculate the short distance
           float distance;
           for(int z = z_lower-bandwise; z <= z_upper+3; z++)
           {
             for(int y = y_lower-3; y < y_upper+3; y++)
             {
               for(int x = x_lower-3; x < x_upper+3; x++)
               {
                 float n = g.gridindex(x, y, z);
                 Vector p = g.pos(x, y, z);
                 Vector p0 = Vector(x1, y1, z1);
                 Vector p1 = Vector(x2, y2, z2);
                 Vector p2 = Vector(x3, y3, z3);

                 Vector temp = p - p0;
                 Vector e0 = p1-p0;
                 Vector e1 = p2-p0;
                 Vector e2 = p2-p1;
                 double e0e1 = e0*e1;
                 double e1e1 = e1*e1;
                 double e0e0 = e0*e0;
                 if(e0e1==0){/*cout<<"divided by zero fault!!index:0"<<endl;*/e0e1=0.00001;}
                 double solution = temp * e0 * e1e1 / e0e1 - temp * e1;
                 double u = (e0e0*e1e1/e0e1 - e0e1);
                 if(u==0){/*cout<<"divided by zero fault!!index:1"<<endl;*/u=0.00001;}
                 u = solution / u;
                 double v = (temp*e0-u*e0e0)/e0e1;
                 if(u>=0&&u<=1&&v>=0&&v<=1&&u+v>=0&&u+v<=1)
                 {
                   Vector d = Vector(x1, y1, z1)+u*e0+v*e1;
                   distance = (p-d).magnitude();
                 }
                 else
                 {
                   if(e0e0==0){/*cout<<"divided by zero fault!!index:2"<<endl;*/e0e0=0.00001;}
                   double umin = -1 * e0 * (p0 - p) / e0e0;
                   if(umin < 0)
                   {
                     umin = 0;
                   }
                   if(umin > 1)
                   {
                     umin = 1;
                   }
                   Vector v01 = (p0 - p + umin * e0);
                   double d01 = v01 * v01;

                   if(e1e1==0){/*cout<<"divided by zero fault!!index:3"<<endl;*/e1e1=0.00001;}
                   umin = -1 * e1 * (p0 - p) / e1e1;
                   if(umin < 0)
                   {
                     umin = 0;
                   }
                   if(umin > 1)
                   {
                     umin = 1;
                   }
                   Vector v02 = (p0 - p + umin * e1);
                   double d02 = v02 * v02;

                   float e2e2 = e2 * e2;
                   if(e2e2==0){/*cout<<"divided by zero fault!!index:4"<<endl;*/e2e2=0.00001;}
                   umin = -1 * e2 * (p1 - p) / e2e2;
                   if(umin < 0)
                   {
                     umin = 0;
                   }
                   if(umin > 1)
                   {
                     umin = 1;
                   }
                   Vector v12 = (p1 - p + umin * e2);
                   double d12 = v12 * v12;
                   double min_dis = min(d01, d02);
                   min_dis = min(min_dis, d12);
                   distance = min_dis;
                 }

                 g.get_data()[n] = min(distance, g.get_data()[n]);
                 //cout<<"gridvalue: "<<g.get_data()[n]<<endl;
               }
             }
           }
//cout<<"test3333333333:  "<<test<<endl;


         }

         //cout<<"test11111111111111:  "<<test<<endl;
         //for all grids, test intersections with triangles
         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               int intersect_num = 0;
               Vector p = g.pos(x, y, z);
               float n = g.gridindex(x, y, z);
               if(g.get_data()[n]<100000)
               {
                 //loop all triangles
                 for(int i=0; i<size; i++)
                 {
                   int point1, point2, point3;
                   point1 = f[i][0]-1;
                   point2 = f[i][1]-1;
                   point3 = f[i][2]-1;
                   double x1, y1, z1;
                   x1 = vp[point1][0];
                   y1 = vp[point1][1];
                   z1 = vp[point1][2];
                   double x2, y2, z2;
                   x2 = vp[point2][0];
                   y2 = vp[point2][1];
                   z2 = vp[point2][2];
                   double x3, y3, z3;
                   x3 = vp[point3][0];
                   y3 = vp[point3][1];
                   z3 = vp[point3][2];

                   Vector p0 = Vector(x1, y1, z1);
                   Vector p1 = Vector(x2, y2, z2);
                   Vector p2 = Vector(x3, y3, z3);
                   //Vector temp = p - p0;
                   Vector e0 = p1-p0;
                   Vector e1 = p2-p0;
                   //Vector e2 = p2-p1;
                   //double e0e1 = e0*e1;
                   //double e1e1 = e1*e1;
                   //double e0e0 = e0*e0;

                   Vector ray_dir(1.0, 0, 0);
                   float r_mult_e0e1 = ray_dir * (e0 ^ e1);
                   if(r_mult_e0e1==0){/*cout<<"divided by zero fault!!index:5"<<endl;*/r_mult_e0e1=0.00001;}
                   float t = (p0 - p)*(e0 ^ e1) / r_mult_e0e1;
                   if(t<0)
                   {

                   }
                   else
                   {
                     Vector e1e0 = e1 ^ e0;
                     Vector e0e1 = e0 ^ e1;
                     Vector p_minus_p0 = p + t*ray_dir - p0;
                     float temp1, temp2;
                     temp1 = e1e0*e1e0;
                     temp2 = e0e1*e0e1;
                     if(temp1==0){/*cout<<"divided by zero fault!!index:6"<<endl;*/temp1=0.00001;}
                     if(temp2==0){/*cout<<"divided by zero fault!!index:7"<<endl;*/temp2=0.00001;}
                     double u = (e1 ^ p_minus_p0) * e1e0 / temp1;
                     double v = (e0 ^ p_minus_p0) * e0e1 / temp2;
                     if(u>=0.0&&u<=1.0&&v>=0.0&&v<=1.0&&u+v<=1.0)
                     {
                       intersect_num++;
                     }
                   }
                 }
                 if(intersect_num%2==0)
                 {
                   g.get_data()[n] = -1 * g.get_data()[n];
                 }
                 // else
                 // {
                 //   std::cout<<"positive"<<std::endl;
                 //   do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                 // }
               }

                //= dsm_ijk(p);
             }
           }
         }

         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               int intersect_num = 0;
               Vector p = g.pos(x, y, z);
               float n = g.gridindex(x, y, z);
               if(g.get_data()[n]<100000)
               {
                 //loop all triangles
                 for(int i=0; i<size; i++)
                 {
                   int point1, point2, point3;
                   point1 = f[i][0]-1;
                   point2 = f[i][1]-1;
                   point3 = f[i][2]-1;
                   double x1, y1, z1;
                   x1 = vp[point1][0];
                   y1 = vp[point1][1];
                   z1 = vp[point1][2];
                   double x2, y2, z2;
                   x2 = vp[point2][0];
                   y2 = vp[point2][1];
                   z2 = vp[point2][2];
                   double x3, y3, z3;
                   x3 = vp[point3][0];
                   y3 = vp[point3][1];
                   z3 = vp[point3][2];

                   Vector p0 = Vector(x1, y1, z1);
                   Vector p1 = Vector(x2, y2, z2);
                   Vector p2 = Vector(x3, y3, z3);
                   Vector e0 = p1-p0;
                   Vector e1 = p2-p0;

                   Vector ray_dir(0.0, 1.0, 0.0);
                   float r_mult_e0e1 = ray_dir * (e0 ^ e1);
                   if(r_mult_e0e1==0){/*cout<<"divided by zero fault!!index:5"<<endl;*/r_mult_e0e1=0.00001;}
                   float t = (p0 - p)*(e0 ^ e1) / r_mult_e0e1;
                   if(t<0)
                   {

                   }
                   else
                   {
                     Vector e1e0 = e1 ^ e0;
                     Vector e0e1 = e0 ^ e1;
                     Vector p_minus_p0 = p + t*ray_dir - p0;
                     float temp1, temp2;
                     temp1 = e1e0*e1e0;
                     temp2 = e0e1*e0e1;
                     if(temp1==0){/*cout<<"divided by zero fault!!index:6"<<endl;*/temp1=0.00001;}
                     if(temp2==0){/*cout<<"divided by zero fault!!index:7"<<endl;*/temp2=0.00001;}
                     double u = (e1 ^ p_minus_p0) * e1e0 / temp1;
                     double v = (e0 ^ p_minus_p0) * e0e1 / temp2;
                     if(u>=0.0&&u<=1.0&&v>=0.0&&v<=1.0&&u+v<=1.0)
                     {
                       intersect_num++;
                     }
                   }
                 }
                 if(intersect_num%2==0)
                 {
                   if(g.get_data()[n]>0)
                      g.get_data()[n] = -1 * g.get_data()[n];
                 }
                 // else
                 // {
                 //   std::cout<<"positive"<<std::endl;
                 //   do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                 // }
               }

             }
           }
         }


         //update grids outside of the band
         //for all grids
         bool empty = true;

         //cout<<"test22222222222222222:  "<<test<<endl;
         while(empty)
         {
           empty = false;
           for(int z=0; z< g.Z(); z++)
           {
             for(int y = 0; y < g.Y(); y++)
             {
               for(int x = 0; x< g.X(); x++)
               {
                 //Vector p = g.pos(x, y, z);
                 float n = g.gridindex(x, y, z);
                 // if(n==0&&g.get_data()[0]!=1000000.0)
                 // {
                 //   std::cout<<"data: "<<g.get_data()[0]<<std::endl;
                 //   do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                 // }


                 //find grids not updated yet
                 if(g.get_data()[n]==1000000.0)//>=1000000)
                 {
                   empty = true;
                   //g.get_data()[n] = 0;//method 1, wrong
                   //method2:
                   //----------------------levelset update(6 directions)---------------------
                   if(x-1>0)
                   {
                     float n1 = g.gridindex(x-1, y, z);
                     if(g.get_data()[n1]!=1000000.0)
                     {
                       if(g.get_data()[n1]>0)
                       {
                         g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                       }
                       else
                       {
                         g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                       }
                       // std::cout<<"data: "<<g.get_data()[n]<<std::endl;
                       // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                       continue;
                     }
                   }

                   if(x+1<g.X())
                   {
                     float n1 = g.gridindex(x+1, y, z);
                     if(g.get_data()[n1]!=1000000.0)
                     {
                       if(g.get_data()[n1]>0)
                       {
                         g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                       }
                       else
                       {
                         g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                       }
                       continue;
                      }
                    }

                    if(y-1>0)
                    {
                      float n1 = g.gridindex(x, y-1, z);
                      if(g.get_data()[n1]!=1000000.0)
                      {
                        if(g.get_data()[n1]>0)
                        {
                          g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                        }
                        else
                        {
                          g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                        }
                        // std::cout<<"data: "<<g.get_data()[n]<<std::endl;
                        // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                        continue;
                       }
                     }

                   if(y+1<g.Y())
                   {
                     float n1 = g.gridindex(x, y+1, z);
                     if(g.get_data()[n1]!=1000000.0)
                     {
                       if(g.get_data()[n1]>0)
                       {
                         g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                       }
                       else
                       {
                         g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                       }
                       continue;
                      }
                    }

                    if(z-1>0)
                    {
                      float n1 = g.gridindex(x, y, z-1);
                      if(g.get_data()[n1]!=1000000.0)
                      {
                        if(g.get_data()[n1]>0)
                        {
                          g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                        }
                        else
                        {
                          g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                        }
                        continue;
                       }
                     }

                   if(z+1<g.Z())
                   {
                     float n1 = g.gridindex(x, y, z+1);
                     if(g.get_data()[n1]!=1000000.0)
                     {
                       if(g.get_data()[n1]>0)
                       {
                         g.get_data()[n] = g.get_data()[n1] + g.get_delta();
                       }
                       else
                       {
                         g.get_data()[n] = g.get_data()[n1] - g.get_delta()*0.01;
                       }
                       continue;
                      }
                    }

                  }
                 // else if(g.get_data()[n]>=0)
                 // {
                 //   // cout<<"data: "<<"n "<<n<<" "<<g.get_data()[n]<<endl;
                 //   // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
                 //   g.get_data()[n] = 1;
                 // }



               }
             }
           }
         }

         for(int z=0; z< g.Z(); z++)
         {
           for(int y = 0; y < g.Y(); y++)
           {
             for(int x = 0; x< g.X(); x++)
             {
               float n = g.gridindex(x, y, z);
               if(g.get_data()[n]>0)
                 g.get_data()[n]=1.0;
               else
                 g.get_data()[n]=0;
               //if(g.get_data()[n]>0)
               //g.get_data()[n]=g.get_data()[n]*10;
               //cout<<"gridvalue: "<<g.get_data()[n]<<endl;

              }
            }
          }

       }

       ~LS_Grids(){}



       Grids get_grid()
       {
         return g;
       }

       private:

       //GridsToFields *gtf;
       Grids g;
       OBJ_model obj;
       int bandwise = 3;
    };
//===================================levelset=========================================


}
#endif

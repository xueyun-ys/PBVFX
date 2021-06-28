

//
//
// Wedge.h
//
// Generate wedges for different noise types
//
//
//


#ifndef __WEDGE_H__
#define __WEDGE_H__

#include "Vector.h"
#include "grids.h"
#include <vector>
#include <algorithm>
#include "Noise.h"
#include "PerlinNoise.h"
#include <iostream>
#include <time.h>
#include <random>

using namespace std;

namespace lux
{

class NoiseWedge
{
  public:

    NoiseWedge(Grids G, Vector O, float r, Noise_t mynoise) :
       g (G),
       radius (r)
    {
      //float s = g.get_delta();
      FractalSum<PerlinNoise> *fs = new FractalSum<PerlinNoise>;

      for(int z=0; z< g.Z(); z++)
      {
          for(int y = 0; y < g.Y(); y++)
          {
            for(int x = 0; x< g.X(); x++)
            {
              float n = g.gridindex(x, y, z);
              g.get_data()[n] = 0.0;
            }
          }
       }

       // struct Noise_t noise_a;
       // noise_a.falloff = 0.8;
       // noise_a.octaves = 3.0;
       // noise_a.frequency = 0.666666;
       // noise_a.axis = Vector(0,0,1);
       // noise_a.amplitude = 3;
       fs->setParameters(mynoise);
       //struct Noise_t noise_b;
       //fs->getParameters(noise_b);
       //cout<<"test: "<<noise_b.frequency<<" "<<noise_b.octaves<<endl;

      for(int z=0; z< g.Z(); z++)
      {
          for(int y = 0; y < g.Y(); y++)
          {
            for(int x = 0; x< g.X(); x++)
            {
              Vector p = g.pos(x, y, z);
              float n = g.gridindex(x, y, z);
              float distance;
              distance = (p - O).magnitude();
              if(distance<r)
              {
                if(distance<mynoise.falloff*r)
                  g.get_data()[n] = max((float)0.0, fs->eval(p));// + 1.0 - distance;
                else
                  g.get_data()[n] = max((float)0.0, fs->eval(p)) * (1-distance/r)/(1-mynoise.falloff);

              }
            }
          }
       }

    }

   ~NoiseWedge(){}

   Grids get_grid()
   {
     return g;
   }



  private:

    Grids g;
    float radius;

};


class PyroclasticWedge
{
  public:

    PyroclasticWedge(Grids G, float R) :
       g (G),
       R (R)
    {
      //float s = g.get_delta();

       noise_a.frequency = 0.666666;//
       noise_a.axis = Vector(0,0,1);
       noise_a.falloff = 0.8;
       noise_a.gamma = 0.8;//
       noise_a.octaves = 3.0;//
       noise_a.amplitude = 2;
       noise_a.fjump = 3;//
       fs->setParameters(noise_a);
    }

   ~PyroclasticWedge(){}

   void Setnoise(struct Noise_t noise_b)//(struct)Noise_t
   {
     noise_a = noise_b;
   }
   void eval()
   {
     for(int z=0; z< g.Z(); z++)
     {
         for(int y = 0; y < g.Y(); y++)
         {
           for(int x = 0; x< g.X(); x++)
           {
             Vector p = g.pos(x, y, z);
             float n = g.gridindex(x, y, z);
             float pn = fs->eval(p.unitvector() * R);
             pn *= 2.0;
             if(pn>1)
                pn = 1;
             //    cout<<"wow--------"<<endl;
             float temp = R - p.magnitude() + pow(abs(pn),  noise_a.gamma);
             if(temp < 0)
               temp = 0;
             else
               temp = 1.0;
             g.get_data()[n] = temp;

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
    float R;
    struct Noise_t noise_a;
    FractalSum<PerlinNoise> *fs = new FractalSum<PerlinNoise>;
};

class WispsWedge
{
  public:

    WispsWedge(Grids G,  float pscale, Vector O, float rho, struct Noise_t noise) ://grids, radius, center of sphere  //, float clump
       g (G),
       pscale (pscale),
       rho(rho),
       noise_a(noise)
       //clump(clump)
    {
      for(int z=0; z< g.Z(); z++)
      {
          for(int y = 0; y < g.Y(); y++)
          {
            for(int x = 0; x< g.X(); x++)
            {
              float n = g.gridindex(x, y, z);
              g.get_data()[n] = 0.0;
            }
          }
       }

       //Get a random radius based on FSPN
       FractalSum<PerlinNoise> *fs1 = new FractalSum<PerlinNoise>;

       struct Noise_t noise_c;
       noise_c.roughness = 0.5;
       noise_c.gamma = 1.5;
       noise_c.octaves = 5.0;//
       noise_c.frequency = 3.0;//
       //noise_a.axis = Vector(0,0,1);
       noise_c.amplitude = 3;
       noise_c.fjump = 2.5;//around 2
       fs1->setParameters(noise_c);

       FractalSum<PerlinNoise> *fs2 = new FractalSum<PerlinNoise>;
       struct Noise_t noise_b;
       //noise_b.falloff = 1.1;
       noise_b.roughness = 0.8;//0-1
       noise_b.gamma = 0.8;
       noise_b.octaves = 5.0;
       noise_b.frequency = 0.5;
       //noise_b.axis = Vector(0,0,1);
       noise_b.amplitude = 1;
       noise_b.fjump = 2;
       fs2->setParameters(noise_b);

      for(int l = 0;l<3000000;l++)
      {
         //srand(time(0));
         std::random_device rd;
         std::mt19937 gen(rd());
         std::uniform_real_distribution<> dis(-1, 1);//uniform distribution between 0 and 1
         // for(int n=0;n<10;n++)
         // {
            // cout<<dis(gen)<<endl;//' '
            // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
         // }
         //pick a random point in unit cube
         float prng_x = dis(gen);
         float prng_y = dis(gen);
         float prng_z = dis(gen);
         // cout<<prng_x<<' '<<prng_y<<' '<<prng_z<<' '<<endl;//' '
         // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

         //project to unit sphere
         float temp = sqrt(prng_x*prng_x + prng_y*prng_y + prng_z*prng_z);
         if(temp==0)
         {
           temp = 0.001;
         }
         float sphere_x = prng_x/temp;
         float sphere_y = prng_y/temp;
         float sphere_z = prng_z/temp;
         Vector prn(prng_x, prng_y, prng_z);

         //Get a random radius based on FSPN
         // FractalSum<PerlinNoise> *fs1 = new FractalSum<PerlinNoise>;
         // struct Noise_t noise_a;
         // noise_a.roughness = 0.5;
         // noise_a.gamma = 1.5;
         // noise_a.octaves = 5.0;//
         // noise_a.frequency = 3.0;//
         // //noise_a.axis = Vector(0,0,1);
         // noise_a.amplitude = 3;
         // noise_a.fjump = 2.5;//around 2
         // fs1->setParameters(noise_a);
         float r = pow( abs(fs1->eval(prn)),  clump);

         //multiply radius
         float xx = r * sphere_x;
         float yy = r * sphere_y;
         float zz = r * sphere_z;
         Vector fspn1(xx, yy, zz);
         Vector fspn2(xx+shift, yy+shift, zz+shift);
         Vector fspn3(xx-shift, yy-shift, zz-shift);


         //Transform to our space
         float xxx = xx * pscale + O.X();
         float yyy = yy * pscale + O.Y();
         float zzz = zz * pscale + O.Z();

         //Vector displacement
         // FractalSum<PerlinNoise> *fs2 = new FractalSum<PerlinNoise>;
         // struct Noise_t noise_b;
         // //noise_b.falloff = 1.1;
         // noise_b.roughness = 0.8;//0-1
         // noise_b.gamma = 0.8;
         // noise_b.octaves = 5.0;
         // noise_b.frequency = 0.5;
         // //noise_b.axis = Vector(0,0,1);
         // noise_b.amplitude = 1;
         // noise_b.fjump = 2;
         // fs2->setParameters(noise_b);
         // // fs2->getParameters(noise_a);
         // // cout<<"roughness: "<<noise_a.roughness<<endl;
         float dx = fs2->eval(fspn1);
         float dy = fs2->eval(fspn2);
         float dz = fs2->eval(fspn3);

         //Use displacement
         float xxxx = xxx + dx;
         float yyyy = yyy + dy;
         float zzzz = zzz + dz;
         //new position
         Vector q(xxxx, yyyy, 1*zzzz);

         //Bake a wisp dot at q in the grid
         Vector llc = g.get_llc();
         float step = g.get_delta();
         int xnum = g.X();
         int ynum = g.Y();
         int znum = g.Z();

         //get relative location
         Vector relative_location = q - llc;
         //if comment this part out and do it in grids function "wispbake", will get some nan point value in one random dimension somehow...
         if(relative_location.X()<=(xnum-1)*step && relative_location.X()>=0
         && relative_location.Y()<=(ynum-1)*step && relative_location.Y()>=0
         && relative_location.Z()>=(znum-1)*step*-1 && relative_location.Z()<=0)
         //&& relative_location.Z()>=llc.Z() - (znum-1)*step && relative_location.Z()<=0)
         {
           //find index position
           //cout<<"hello!"<<endl;
           g.BakeWisp(q, rho);

         }
         else
         {
           continue;
         }
      }

      //mask
      // float maxn = -10000;
      // float minn = 10000;
      // for(int z=0; z< g.Z(); z++)
      // {
      //     for(int y = 0; y < g.Y(); y++)
      //     {
      //       for(int x = 0; x< g.X(); x++)
      //       {
      //         //Vector p = g.pos(x, y, z);
      //         float n = g.gridindex(x, y, z);
      //         if(maxn < g.get_data()[n])
      //         {
      //           maxn = g.get_data()[n];
      //         }
      //         if(minn > g.get_data()[n])
      //         {
      //           minn = g.get_data()[n];
      //         }
      //         // if(g.get_data()[n] < 0)
      //         //   cout<<"Wrong! negetive value for wispwedge"<<endl;
      //         //   //g.get_data()[n] = 0;
      //         // else if(g.get_data()[n]>0.8)
      //         //
      //         //   g.get_data()[n] = 1.0;
      //         //   //g.get_data()[n] *= 1000.0;
      //       }
      //     }
      //  }
      //  cout<<"max: "<<maxn<<"min: "<<minn<<endl;
    }

    WispsWedge(Grids G,  float pscale, Vector O, float rho) ://grids, radius, center of sphere
       g (G),
       pscale (pscale),
       rho(rho)
    {
      for(int z=0; z< g.Z(); z++)
      {
          for(int y = 0; y < g.Y(); y++)
          {
            for(int x = 0; x< g.X(); x++)
            {
              float n = g.gridindex(x, y, z);
              g.get_data()[n] = 0.0;
            }
          }
       }

       //Get a random radius based on FSPN
       FractalSum<PerlinNoise> *fs1 = new FractalSum<PerlinNoise>;

       noise_a.roughness = 0.5;
       noise_a.gamma = 1.5;
       noise_a.octaves = 5.0;//
       noise_a.frequency = 3.0;//
       //noise_a.axis = Vector(0,0,1);
       noise_a.amplitude = 3;
       noise_a.fjump = 2.5;//around 2
       fs1->setParameters(noise_a);

       FractalSum<PerlinNoise> *fs2 = new FractalSum<PerlinNoise>;
       struct Noise_t noise_b;
       //noise_b.falloff = 1.1;
       noise_b.roughness = 0.8;//0-1
       noise_b.gamma = 0.8;
       noise_b.octaves = 5.0;
       noise_b.frequency = 0.5;
       //noise_b.axis = Vector(0,0,1);
       noise_b.amplitude = 1;
       noise_b.fjump = 2;
       fs2->setParameters(noise_b);

      for(int l = 0;l<3000000;l++)
      {
         //srand(time(0));
         std::random_device rd;
         std::mt19937 gen(rd());
         std::uniform_real_distribution<> dis(-1, 1);//uniform distribution between 0 and 1
         //pick a random point in unit cube
         float prng_x = dis(gen);
         float prng_y = dis(gen);
         float prng_z = dis(gen);
         // cout<<prng_x<<' '<<prng_y<<' '<<prng_z<<' '<<endl;//' '
         // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');

         //project to unit sphere
         float temp = sqrt(prng_x*prng_x + prng_y*prng_y + prng_z*prng_z);
         if(temp==0)
         {
           temp = 0.001;
         }
         float sphere_x = prng_x/temp;
         float sphere_y = prng_y/temp;
         float sphere_z = prng_z/temp;
         Vector prn(prng_x, prng_y, prng_z);

         //Get a random radius based on FSPN
         float r = pow( abs(fs1->eval(prn)),  clump);

         //multiply radius
         float xx = r * sphere_x;
         float yy = r * sphere_y;
         float zz = r * sphere_z;
         Vector fspn1(xx, yy, zz);
         Vector fspn2(xx+shift, yy+shift, zz+shift);
         Vector fspn3(xx-shift, yy-shift, zz-shift);


         //Transform to our space
         float xxx = xx * pscale + O.X();
         float yyy = yy * pscale + O.Y();
         float zzz = zz * pscale + O.Z();

         //Vector displacement
         float dx = fs2->eval(fspn1);
         float dy = fs2->eval(fspn2);
         float dz = fs2->eval(fspn3);

         //Use displacement
         float xxxx = xxx + dx;
         float yyyy = yyy + dy;
         float zzzz = zzz + dz;
         //new position
         Vector q(xxxx, yyyy, 1*zzzz);

         //Bake a wisp dot at q in the grid
         Vector llc = g.get_llc();
         float step = g.get_delta();
         int xnum = g.X();
         int ynum = g.Y();
         int znum = g.Z();

         //get relative location
         Vector relative_location = q - llc;
         //if comment this part out and do it in grids function "wispbake", will get some nan point value in one random dimension somehow...
         if(relative_location.X()<=(xnum-1)*step && relative_location.X()>=0
         && relative_location.Y()<=(ynum-1)*step && relative_location.Y()>=0
         && relative_location.Z()>=(znum-1)*step*-1 && relative_location.Z()<=0)
         {

           g.BakeWisp(q, rho);

         }
         else
         {
           continue;
         }
      }

      //mask
      // float maxn = -10000;
      // float minn = 10000;
      // for(int z=0; z< g.Z(); z++)
      // {
      //     for(int y = 0; y < g.Y(); y++)
      //     {
      //       for(int x = 0; x< g.X(); x++)
      //       {
      //         //Vector p = g.pos(x, y, z);
      //         float n = g.gridindex(x, y, z);
      //         if(maxn < g.get_data()[n])
      //         {
      //           maxn = g.get_data()[n];
      //         }
      //         if(minn > g.get_data()[n])
      //         {
      //           minn = g.get_data()[n];
      //         }
      //
      //       }
      //     }
      //  }
    }

   ~WispsWedge(){}

   Grids get_grid()
   {
     return g;
   }



  private:

    Grids g;
    float pscale;
    float clump = 0.3;//0.6;//0.3333
    float shift = 0.12;//around0.1
    float rho = 0.0001;
    struct Noise_t noise_a;
};


}

#endif

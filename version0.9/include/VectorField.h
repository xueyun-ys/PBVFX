

//
//
// VectorField.h
//
// Vector Functions
//
//
//


#ifndef __VECTORFIELD_H__
#define __VECTORFIELD_H__

#include "Vector.h"
#include "Matrix.h"
//#include "grids.h"
#include <vector>
#include <algorithm>
//#include "Noise.h"
//#include "PerlinNoise.h"
#include "Wedge.h"
#include <iostream>
#include <time.h>
#include <random>
#include "Volume.h"
#include <cmath>

using namespace std;


namespace lux
{

    class Identity : public Volume<Vector>
    {
      public:

        Identity()
        {}

       ~Identity(){}

       const Vector eval(const Vector& x)const
       {
         return x;
       }



      private:
    };


    class InverseVector : public Volume<Vector>
    {
      public:

        InverseVector(Volume<Vector> *v) : vec(v) {}

       ~InverseVector(){}

       const Vector eval(const Vector& x)const
       {
         return -1.0*vec->eval(x);
       }
       const Matrix grad(const Vector& x)const
       {
         return -1.0*vec->grad(x);
       }

      private:
        Volume<Vector> *vec;
    };

    class AddVector : public Volume<Vector>
    {
      public:

        AddVector(Volume<Vector> *v1, Volume<Vector> *v2) : vec1(v1), vec2(v2) {}

       ~AddVector(){}

       const Vector eval(const Vector& x)const
       {
         return vec1->eval(x) + vec2->eval(x);
       }
       const Matrix grad(const Vector& x)const
       {
         return vec1->grad(x) + vec2->grad(x);
       }

      private:
        Volume<Vector> *vec1;
        Volume<Vector> *vec2;
    };


    class MinusVector : public Volume<Vector>
    {
      public:

        MinusVector(Volume<Vector> *v1, Volume<Vector> *v2) : vec1(v1), vec2(v2) {}

       ~MinusVector(){}

       const Vector eval(const Vector& x)const
       {
         return vec1->eval(x) - vec2->eval(x);
       }
       const Matrix grad(const Vector& x)const
       {
         return vec1->grad(x) - vec2->grad(x);
       }

      private:
        Volume<Vector> *vec1;
        Volume<Vector> *vec2;
    };

    class InnerproductVector : public Volume<float>
    {
      public:

        InnerproductVector(Volume<Vector> *v1, Volume<Vector> *v2) : vec1(v1), vec2(v2) {}

       ~InnerproductVector(){}

       const float eval(const Vector& x)const
       {
         return vec1->eval(x) * vec2->eval(x);
       }
       // const Vector grad(const Vector& x)const
       // {
       //   return vec1->grad(x) * vec2->eval(x) + vec1->eval(x) * vec2->grad(x);//Matrix * Vector
       // }

      private:
        Volume<Vector> *vec1;
        Volume<Vector> *vec2;
    };
    //class OutterproductVector : public Volume<Vector>

//Vector Scalar Multiple
    class VSMVector : public Volume<Vector>
    {
      public:

        VSMVector(Volume<Vector> *v, Volume<float> *f) : vec(v), sca(f) {}

       ~VSMVector(){}

       const Vector eval(const Vector& x)const
       {
         return vec->eval(x) * sca->eval(x);
       }
       // const Matrix grad(const Vector& x)const
       // {
       //   cout<<"wrong"<<endl;
       //   return vec->grad(x) * sca->eval(x) + vec->eval(x) @outerproduct sca->grad(x);
       // }

      private:
        Volume<Vector> *vec;
        Volume<float> *sca;
    };

    //First Handy-Dandy VF(Vector Field)

    //Second Handy-Dandy VF(Vector Field)
    class GradSF : public Volume<Vector>
    {
      public:

        GradSF(Volume<float> *f) : sca(f) {}

       ~GradSF(){}

       const Vector eval(const Vector& x)const
       {
         return sca->grad(x);
       }
       // const Matrix grad(const Vector& x)const
       // {

       // }
      private:
        Volume<float> *sca;
    };

//=======================================operations==================================

  //Closet Point Transform
  class CPT : public Volume<Vector>
  {
    public:

      CPT(Volume<float> *f) : sca(f) {}

      ~CPT(){}

      const Vector eval(const Vector& x)const
      {
        //return Vector(x.X(), -3.0, x.Z());
        Vector temp = x - sca->eval(x) * sca->grad(x);
        // if(temp.Y()<0.0001)
        // temp=Vector(temp.X(),0.0,temp.Z());
        return temp;
      }

     private:
       Volume<float> *sca;
   };

   //Warp
   class Warp : public Volume<float>
   {
     public:

       Warp(Volume<float> *f, Volume<Vector> *v) : sca(f), vec(v) {}

       ~Warp(){}

       const float eval(const Vector& x)const
       {
         Vector warpx = vec->eval(x);
         return sca->eval(warpx);
       }

      private:
        Volume<Vector> *vec;
        Volume<float> *sca;
    };

//FSPN
    class FSPN : public Volume<float>
    {
      public:

        FSPN(Noise_t N) : noise_a(N) {fs.setParameters(noise_a);}
        FSPN() {fs.setParameters(noise_a);}
        ~FSPN(){}

        const float eval(const Vector& x)const
        {

           //fs->setParameters(noise_a);

           // if(x.Y()<0.0001)
           // x=Vector(x.X(),0.0,x.Z());
           //x.Y()=0;
           Vector X(1.0, 0.0, -1.0);
           float raw = fs.eval(x) * scale;//X
           // if(raw>0)
           // {
           //
           // }
           // else
           // {
           //
           // }
           // cout<<"fspn "<<X.X()<<" "<<X.Y()<<" "<<X.Z()<<endl;
           // cout<<"raw: "<<raw<<endl;
           // do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
           return raw;
        }

       private:
         struct Noise_t noise_a;
         struct Noise_t noise_b;
         float scale = 1.0;
         //FractalSum<PerlinNoiseGustavson> *fs = new FractalSum<PerlinNoiseGustavson>;
         //FractalSum<BaseNoise> *fs = new FractalSum<BaseNoise>;
         FractalSum<PerlinNoise> fs;
     };

     //NPT
     class NPT : public Volume<Vector>
     {
       public:

         NPT(Volume<float> *f) : f(f) {}
         ~NPT(){}

         const Vector eval(const Vector& x)const
         {
            Vector temp = f->grad(x);
            return x - f->eval(x)*temp/(temp*temp);
         }

        private:
          Volume<float> *f;
      };
      class nNPT : public Volume<Vector>
      {
        public:

          nNPT(NPT *npt, int num) : npt(npt), num(num) {}
          ~nNPT(){}

          const Vector eval(const Vector& x)const
          {
            Vector temp = x;
            //Vector p = x;
            for(int n =1;n<num+1;n++)
            {
              temp = npt->eval(temp);

            }
             return temp;
          }

         private:
           Volume<Vector> *npt;
           int num = 0;
       };

      //Pyroclastic Displacement
      class PyroDisplacement : public Volume<float>
      {
        public:

          PyroDisplacement(Volume<float> *f, Volume<float> *warp) : f(f),warp(warp) {}
          ~PyroDisplacement(){}

          const float eval(const Vector& x)const
          {
            //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
              float raw = warp->eval(x);
              if(raw<0)
              {
                float result = -1.0*pow((raw*-1.0), 1.0)*3;
                //cout<<"result: "<<result<<endl;
                return f->eval(x)+result;
              }
              else
              {
                float result = pow(raw, 1.0)*10;
                //cout<<"result: "<<result<<endl;
                return f->eval(x)+result;
              }
             //return f->eval(x) + pow(abs(warp->eval(x)), gamma);
          }

         private:
           Volume<float> *f;
           Volume<float> *warp;
           float gamma = 2.0;
       };


//=======================================operations==================================
}

#endif

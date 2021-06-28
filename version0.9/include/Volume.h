
#ifndef __VOLUME_H__
#define __VOLUME_H__

#include "Vector.h"
#include "Matrix.h"
//#include "Forms.h"
#include "Color.h"
#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include <cstring>

//#include "VolumeGrid.h"
//#include "SparseGrid.h"

namespace lux
{


//-----------------------------------------------------------------------------
// Setting up logic to be able to determine the data type of the gradient
template <typename U>
struct GradType
{
   typedef int GType;
};

template<>
struct GradType<float>
{
   typedef Vector GType;
};

template<>
struct GradType<Vector>
{
   typedef Matrix GType;
};



//-----------------------------------------------------------------------------

template< typename U >
class Volume
{
  public:

    Volume(){}

   virtual ~Volume(){}

   typedef U volumeDataType;
   typedef typename GradType<U>::GType volumeGradType;

   virtual const volumeDataType eval( const Vector& P ) const { volumeDataType base; std::cout<<"ooops! wrong coding1"<<std::endl; return base; }
   virtual const volumeGradType grad( const Vector& P ) const { volumeGradType base; std::cout<<"ooops! wrong coding2"<<std::endl; return base; }
   // virtual const volumeGradType grad( const Vector& P ) const
   // {
   //   //finite difference approach
   //   float delta = 0.01;
   //   volumeDataType fxp = eval(P+Vector(delta, 0, 0));
   //   volumeDataType fxm = eval(P-Vector(delta, 0, 0));
   //   volumeDataType gx = (fxp-fxm)/(2*delta);
   //   volumeDataType fyp = eval(P+Vector(0, delta, 0));
   //   volumeDataType fym = eval(P-Vector(0, delta, 0));
   //   volumeDataType gy = (fyp-fym)/(2*delta);
   //   volumeDataType fzp = eval(P+Vector(0, 0, delta));
   //   volumeDataType fzm = eval(P-Vector(0, 0, delta));
   //   volumeDataType gz = (fzp-fzm)/(2*delta);
   //   //volumeGradType base( gx,  gy,  gz);
   //   return volumeGradType( gx,  gy,  gz);
   //   //return base;
   // }
   //virtual volumeDataType eval( Vector& P ) { volumeDataType base;  return base; }//std::cout<<"ooops! wrong coding3"<<std::endl;
   //virtual volumeGradType grad( Vector& P ) { volumeGradType base; std::cout<<"ooops! wrong coding2"<<std::endl; return base; }


};

typedef Volume<float>* VolumeFloatPtr;
typedef Volume<Color>* VolumeColorPtr;
typedef Volume<Vector>* VolumeVectorPtr;
typedef Volume<Matrix>* VolumeMatrixPtr;



}

#endif

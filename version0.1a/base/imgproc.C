
#include <cmath>
#include "imgproc.h"
#include <string>
#include <iostream>

#include "Objects.h"
#include "Camera.h"
#include "Color.h"
#include "Image.h"

#include <OpenImageIO/imageio.h>
OIIO_NAMESPACE_USING

using namespace img;
using namespace std;
using namespace lux;



ImgProc::ImgProc() :
  Nx (0),
  Ny (0),
  Nc (0),
  Nsize (0),
  single_Nsize (0),
  img_data (nullptr)
{}

ImgProc::~ImgProc()
{
   clear();
}

void ImgProc::clear()
{
   if( img_data != nullptr ){ delete[] img_data; img_data = nullptr;}
   Nx = 0;
   Ny = 0;
   Nc = 0;
   Nsize = 0;
   single_Nsize = 0;
}

void ImgProc::clear(int nX, int nY, int nC)
{
   clear();
   Nx = nX;
   Ny = nY;
   Nc = nC;
   Nsize = (long)Nx * (long)Ny * (long)Nc;
   single_Nsize = (long)Nx * (long)Ny;
   img_data = new float[Nsize];
#pragma omp parallel for
   for(long i=0;i<Nsize;i++){ img_data[i] = 0.0; }
}

bool ImgProc::load( const std::string& filename )
{
   auto in = ImageInput::create (filename);
   if (!in) {return false;}
   ImageSpec spec;
   in->open (filename, spec);
   clear();
   Nx = spec.width;
   Ny = spec.height;
   Nc = spec.nchannels;
   Nsize = (long)Nx * (long)Ny * (long)Nc;
   single_Nsize = (long)Nx * (long)Ny;
   img_data = new float[Nsize];
   in->read_image(TypeDesc::FLOAT, img_data);
   in->close ();
   return true;
}


void ImgProc::value( int i, int j, std::vector<float>& pixel) const
{
   pixel.clear();
   if( img_data == nullptr ){ return; }
   if( i<0 || i>=Nx ){ return; }
   if( j<0 || j>=Ny ){ return; }
   pixel.resize(Nc);
   for( int c=0;c<Nc;c++ )
   {
      pixel[c] = img_data[index(i,j,c)];
   }
   return;
}

void ImgProc::set_value( int i, int j, const std::vector<float>& pixel)
{
   if( img_data == nullptr ){ return; }
   if( i<0 || i>=Nx ){ return; }
   if( j<0 || j>=Ny ){ return; }
   if( Nc > (int)pixel.size() ){ return; }
#pragma omp parallel for
   for( int c=0;c<Nc;c++ )
   {
      img_data[index(i,j,c)] = pixel[c];
   }
   return;
}


ImgProc::ImgProc(const ImgProc& v) :
  Nx (v.Nx),
  Ny (v.Ny),
  Nc (v.Nc),
  Nsize (v.Nsize),
  single_Nsize (v.single_Nsize)
{
   img_data = new float[Nsize];
#pragma omp parallel for
   for( long i=0;i<Nsize;i++){ img_data[i] = v.img_data[i]; }
}

ImgProc& ImgProc::operator=(const ImgProc& v)
{
   if( this == &v ){ return *this; }
   if( Nx != v.Nx || Ny != v.Ny || Nc != v.Nc )
   {
      clear();
      Nx = v.Nx;
      Ny = v.Ny;
      Nc = v.Nc;
      Nsize = v.Nsize;
      single_Nsize = v.single_Nsize;
   }
   img_data = new float[Nsize];
#pragma omp parallel for
   for( long i=0;i<Nsize;i++){ img_data[i] = v.img_data[i]; }
   return *this;
}


void ImgProc::operator*=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] *= v; }
}

void ImgProc::operator/=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] /= v; }
}

void ImgProc::operator+=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] += v; }
}

void ImgProc::operator-=(float v)
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] -= v; }
}

//==================================================================================
//volume rendering
void ImgProc::render()//input: density field, color field
{
   if( img_data == nullptr ){ return; }
   //clear(Nx, Ny, 4);
   //clear(1920, 1080, 4);
   clear(360, 270, 4);
   //cout<<"Nx:    "<<Nx<<endl;
   //cout<<"Ny:    "<<Ny<<endl;
//#pragma omp parallel for
    //render avator

  	// using namespace std;
  	// using namespace lux;

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
    // ObjectCone *d4 = new ObjectCone(6.0, 0.2, Vector(0.0, -1.0, 0.0));
    // ObjectSFRotation *s8 = new ObjectSFRotation(d4, Vector(0.0, 0.0, 1.0), 3.1415*2.0/3.0);
    // ObjectSFTranslate *t1 = new ObjectSFTranslate(s8, Vector(4.0, 1.7, 0.0));
    //
    // ObjectBox *d5 = new ObjectBox(1.0, 5.0);
    // ObjectSFRotation *t2 = new ObjectSFRotation(d5, Vector(0.0, 0.0, 1.0), 3.1415/6.0);
    // ObjectSFTranslate *t3 = new ObjectSFTranslate(t2, Vector(5.0, 2.5, 0.0));
    // ObjectCombined *s9 = new ObjectCombined(t3, t1);
    //
    // ObjectCone *s10 = new ObjectCone(3.0, 0.4, Vector(0.0, -1.0, 0.0));
    // ObjectSFRotation *s11 = new ObjectSFRotation(s10, Vector(0.0, 0.0, 1.0), -3.1415/3.0);
    // ObjectSFScale *s12 = new ObjectSFScale(s11, Vector(0.0, 0.0, 0.0), 0.7);
    // ObjectSFTranslate *t4 = new ObjectSFTranslate(s12, Vector(5.0, 2.5, 0.0));
    //
    // ObjectCombined *s13 = new ObjectCombined(s9, t4);
    // ObjectSFScale *s14 = new ObjectSFScale(s13, Vector(0.0, 0.0, 0.0), 0.4);
    //
    // ObjectTorus *s15 = new ObjectTorus(Vector(0.0, 0.0, 1.0), 5.0, 1.0);
    // ObjectPlane *s16 = new ObjectPlane(Vector(0.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0));
    // ObjectCutout *s17 = new ObjectCutout(s15, s16);
    // ObjectPlane *s18 = new ObjectPlane(Vector(0.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0));
    // ObjectSFTranslate *s19 = new ObjectSFTranslate(s18, Vector(0.0, -2.0, 0.0));
    // ObjectCutout *s20 = new ObjectCutout(s17, s19);
    // ObjectSFScale *s21 = new ObjectSFScale(s20, Vector(0.0, 0.0, 0.0), 0.27);
    // ObjectSFRotation *t5 = new ObjectSFRotation(s21, Vector(0.0, 0.0, 1.0), 3.1415/6.0);
    //
    // ObjectSFTranslate *s22 = new ObjectSFTranslate(t5, Vector(4.2, 0.6, -0.2));
    // ObjectSFTranslate *s23 = new ObjectSFTranslate(t5, Vector(4.2, 0.6, 0.2));
    // ObjectCombined *s24 = new ObjectCombined(s22, s23);
    //
    // ObjectCombined *s26 = new ObjectCombined(s14, s24);
    // ObjectSFTranslate *s27 = new ObjectSFTranslate(s26, Vector(0.0, 4.0, 0.0));
    // ObjectSFRotation *s28 = new ObjectSFRotation(s27, Vector(0.0, 1.0, 0.0), 3.1415/2.0);
    // ObjectMask *m4 = new ObjectMask(s28);
    //
    // Color c3(0.2d,0.6d,0.3d,0.0d);
    // ConstantColor *colorVolume3 = new ConstantColor(c3);
    // ObjectColorDensityFieldMultiple *cf3 = new ObjectColorDensityFieldMultiple(m4, colorVolume3);
    //
    // ObjectColorFieldAdd *color1 = new ObjectColorFieldAdd(cf3, s2);
    //
    // ObjectDensityFieldAdd *m5 = new ObjectDensityFieldAdd(m1, m4);


    //body
    // ObjectBox *s29 = new ObjectBox(1.0, 5.0);
    // ObjectSFScale2 *s30 = new ObjectSFScale2(s29, 0.7, 2);
    // ObjectSFScale *s31 = new ObjectSFScale(s30, Vector(0.0, 0.0, 0.0), 3.0);
    // ObjectSFTranslate *s32 = new ObjectSFTranslate(s31, Vector(1.0, -1.3, 1.0));
    // //ObjectCutout *s33 = new ObjectCutout(s32, s1);
    //
    // ObjectMask *m6 = new ObjectMask(s32);
    // Color c4(0.4d,0.6d,0.7d,0.0d);
    // ConstantColor *colorVolume4 = new ConstantColor(c4);
    // ObjectColorDensityFieldMultiple *cf4 = new ObjectColorDensityFieldMultiple(m6, colorVolume4);
    //
    // ObjectColorFieldAdd *color2 = new ObjectColorFieldAdd(cf4, s2);//color1);
    // ObjectDensityFieldAdd *m7 = new ObjectDensityFieldAdd(m1, m6);//m5

//feet
    // ObjectIcosahedron *s34 = new ObjectIcosahedron(1.0);
    // ObjectSFScale *s35 = new ObjectSFScale(s34, Vector(0.0, 0.0, 0.0), 0.2);
    // ObjectSFTranslate *s36 = new ObjectSFTranslate(s35, Vector(0.0, -5.0, 0.0));
    // ObjectSFTranslate *s37 = new ObjectSFTranslate(s35, Vector(2.8, -5.0, 0.0));
    // ObjectCombined *s38 = new ObjectCombined(s36, s37);
    // //ObjectCutout *s39 = new ObjectCutout(s38, s32);
    //
    // ObjectMask *m8 = new ObjectMask(s38);
    // Color c5(0.4d,0.8d,0.7d,0.0d);
    // ConstantColor *colorVolume5 = new ConstantColor(c5);
    // ObjectColorDensityFieldMultiple *cf5 = new ObjectColorDensityFieldMultiple(m8, colorVolume5);
    //
    // ObjectColorFieldAdd *color3 = new ObjectColorFieldAdd(cf5, color2);
    // ObjectDensityFieldAdd *m9 = new ObjectDensityFieldAdd(m7, m8);

//arms
    ObjectEllipse *s40 = new ObjectEllipse(Vector(1.0, 0.0, 0.0), 5.0, 2.0 );
    ObjectSFScale *s41 = new ObjectSFScale(s40, Vector(0.0, 0.0, 0.0), 0.5);
    ObjectSFRotation *ts44 = new ObjectSFRotation(s41, Vector(0.0, 0.0, 1.0), -3.1415/3.0);
    ObjectSFRotation *s43 = new ObjectSFRotation(s41, Vector(0.0, 0.0, 1.0), 3.1415/3.0);
    ObjectSFTranslate *s44 = new ObjectSFTranslate(s43, Vector(5.3, 0.0, 0.0));
    ObjectSFTranslate *s42 = new ObjectSFTranslate(ts44, Vector(-3.5, 0.0, 0.0));

    //ObjectSFTranslate *s37 = new ObjectSFTranslate(s35, Vector(2.8, -5.0, 0.0));
    ObjectCombined *s45 = new ObjectCombined(s42, s44);
    //ObjectCutout *s39 = new ObjectCutout(s38, s32);

    ObjectMask *m10 = new ObjectMask(s45);
    Color c6(0.6d,0.5d,0.7d,0.0d);
    ConstantColor *colorVolume6 = new ConstantColor(c6);
    ObjectColorDensityFieldMultiple *cf6 = new ObjectColorDensityFieldMultiple(m10, colorVolume6);

    ObjectColorFieldAdd *color4 = new ObjectColorFieldAdd(cf6, s2);//color3
    ObjectDensityFieldAdd *m11 = new ObjectDensityFieldAdd(m1, m10);//m9

    //ObjectPlane *densityVolume = new ObjectPlane(Vector(0.0, 1.0, 0.0), Vector(1.0, 0.0, 0.0));//????
//weapon
    ObjectCylinder *s46 = new ObjectCylinder(0.7, Vector(0.0, 1.0, 0.0));
    ObjectSFTranslate *s47 = new ObjectSFTranslate(s46, Vector(-4.0, 0.0, 0.0));
    //ObjectSphere *s48 = new ObjectSphere(20.0);
    //ObjectIntersection *s49 = new ObjectIntersection(s47, s48);
    ObjectPlane *s48 = new ObjectPlane(Vector(0.0, 60.0, 0.0), Vector(0.0, -1.0, 0.0));
    ObjectPlane *s49 = new ObjectPlane(Vector(0.0, -60.0, 0.0), Vector(0.0, 1.0, 0.0));
    ObjectCutout *s50 = new ObjectCutout(s47, s48);
    ObjectCutout *s51 = new ObjectCutout(s50, s49);

    ObjectMask *m12 = new ObjectMask(s51);
    Color c7(0.3d,0.2d,0.2d,0.0d);
    ConstantColor *colorVolume7 = new ConstantColor(c7);
    ObjectColorDensityFieldMultiple *cf7 = new ObjectColorDensityFieldMultiple(m12, colorVolume7);

    ObjectColorFieldAdd *color5 = new ObjectColorFieldAdd(cf7, color4);
    ObjectDensityFieldAdd *m13 = new ObjectDensityFieldAdd(m11, m12);


    ObjectSteinerPatch *s52 = new ObjectSteinerPatch();
    ObjectSFTranslate *s53 = new ObjectSFTranslate(s52, Vector(-0.0, 0.0, 0.0));

    ObjectMask *m14 = new ObjectMask(s53);
    Color c8(0.5d,0.3d,0.7d,0.0d);
    ConstantColor *colorVolume8 = new ConstantColor(c8);
    ObjectColorDensityFieldMultiple *cf8 = new ObjectColorDensityFieldMultiple(m14, colorVolume8);

    ObjectColorFieldAdd *color6 = new ObjectColorFieldAdd(cf8, color5);
    ObjectDensityFieldAdd *m15 = new ObjectDensityFieldAdd(m13, m14);

    //ObjectShell *s4 = new ObjectShell(densityVolume2, 0.2);


    int test = 0;

  	Camera cam;
  	//Image img;
    #pragma omp parallel for
  	//outloop should be Ny!
    for (int v = 0; v< Ny ;v++)
  	{
      #pragma omp parallel for
  		for (int u = 0;u< Nx;u++)
  		{
  			Vector Np;
  			float T = 1.0;
  			float K = 0.9;
  			double del_s = 0.8;//0.01
  			double s;
  			float alpha;
  			float del_t;
  			//cam.setEyeViewUp(Vector(0, -0, -2.0), Vector(0, 0, 1.0), Vector(0, 0, -1));
  			//cam.setRight();
  			//cam.setAspectRatio(10);
  			//cam.setFov(178);
        double x = (double)u / (double)Nx;
        double y = (double)v / (double)Ny;
  			Np = cam.view(x, y);
  			//Vector Xp = cam.eye() + u * cam.right() + v * cam.up();
  			//Vector Np = (Xp - cam.eye()) / (cam.view() * (Xp - cam.eye())) - cam.view();
  			Vector X = cam.eye();// + Np * cam.nearPlane();
  			s = cam.nearPlane();
        //cout<<"near:    "<<s<<endl;
        //cout<<"channel#:      "<<Nc<<endl;
  			 //Lp = (0, 0, 0, 0)
  			Color Lp;
        // double temp =densityVolume->eval(X);
        // long ind = index(u,v,0);
        // cout<<"pixel:  "<<ind<<":    "<<endl;
        // cout<<"density:"<<": "<<temp<<endl;
        // double temp2 =s1->eval(X);
        // cout<<"s1:  "<<":    "<<temp2<<endl;
        // double temp3 =s2->eval(X).Z();
        // cout<<"s2:  "<<":    "<<temp3<<endl;

  			while (s < cam.farPlane() && T > 0.0001)
  			{
          //cout<<"do it once: ============================="<<endl;
  				X = X + del_s * Np;
          // if(u==500&&v==360)
          // {
            // cout<<"X!!!!!======================    "<<X.X()<<endl;
            // cout<<"Y!!!!!======================    "<<X.Y()<<endl;
            // cout<<"Z!!!!!======================    "<<X.Z()<<endl;
            // cout<<"s1!!!!!======================    "<<s1->eval(X)<<endl; //1
            // cout<<"s2!!!!!======================    "<<s2->eval(X).Z()<<endl; //0
            // cout<<"colorVolume!!!!!======================    "<<colorVolume->eval(X).Z()<<endl; //0.5
          // }
  				s += del_s;
  				del_t = exp(-1.0 * K * del_s * m13->eval(X));
          //cout<<"1111111111111111111:    "<<endl;
          //double test = s2->eval(X).Z();
  				Lp = Lp + color5->eval(X) * T * (1.0 - del_t) / K;
          //cout<<"2222222222222222222:    "<<endl;
          //system("Pause");
          // if(s2->eval(X).Z()>0)
          // {
          //   cout<<"Cool!!!!!======================    "<<endl;
          // }
          // if(u==500&&v==360)
          // {
          //   cout<<"Lp1:    "<<Lp.Z()<<endl;
          //   cout<<"eval:    "<< test <<endl;
          //   //cout<<"!!!:    "<<Lp.Z()<<endl;
          // }
  				T *= del_t;
          alpha = 1.0 - T;
  				Lp.set(Lp.X(), Lp.Y(), Lp.Z(), alpha);
          //if(Lp.Z()>0)
            //cout<<"!!!:    "<<Lp.Z()<<endl;
  			}
        long dex = index(u,v,0);
        //img_data[index(i,j,c)] = pixel[c];
        img_data[dex] = Lp.X();
        img_data[dex+1] = Lp.Y();
        img_data[dex+2] = Lp.Z();
        img_data[dex+3] = Lp.W();//alpha;
        test ++;
        //if(Lp.Z()>0)
          //cout<<"!!!:    "<<Lp.Z()<<endl;
        cout<<"test:    "<<test<<endl;
  		}

    }
}
//==================================================================================

void ImgProc::compliment()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ ){ img_data[i] = 1.0 - img_data[i]; }


//====undo=====
   undo_vec.push_back(10);
}

//======================================================================================
void ImgProc::brightness_up()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = img_data[i] * 1.1;
   }

//====undo=====
   undo_vec.push_back(1);
}

void ImgProc::brightness_down()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = img_data[i] / 1.1;
   }

//====undo=====
   undo_vec.push_back(2);
}

void ImgProc::bias_up()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = img_data[i] + 0.05;
   }

//====undo=====
   undo_vec.push_back(3);
}

void ImgProc::bias_down()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = img_data[i] - 0.05;
   }

//====undo=====
   undo_vec.push_back(4);
}

void ImgProc::gamma_up()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = pow(img_data[i], 0.9);
   }

//====undo=====
   undo_vec.push_back(5);
}

void ImgProc::gamma_down()
{
   if( img_data == nullptr ){ return; }
#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
       img_data[i] = pow(img_data[i], 1.1);
   }

//====undo=====
   undo_vec.push_back(6);
}

void ImgProc::grayscale()
{
   if( img_data == nullptr ){ return; }
   vector<float> undo_temp;
   undo_temp.resize(Nsize);

#pragma omp parallel for
   for( long i=0;i<Nsize;i = i+Nc )
   {
       for(int j=0; j<Nc; j++)
       {
         undo_temp[i+j] = img_data[i+j];
       }
      float g = img_data[i] * 0.2126 + img_data[i+1] * 0.7152 + img_data[i+2] * 0.0722;
       img_data[i] =  g;
       img_data[i+1] =  g;
       img_data[i+2] =  g;
   }

//====undo=====
   undo_vec.push_back(7);
   temp_data.push_back(undo_temp);
}
//pixel[c] = img_data[index(i,j,c)];

void ImgProc::quantize()
{
   if( img_data == nullptr ){ return; }
   vector<float> undo_temp;
   undo_temp.resize(Nsize);

#pragma omp parallel for
   for( long i=0;i<Nsize;i++ )
   {
      undo_temp[i] = img_data[i];
      img_data[i] = floor(img_data[i] * 10) / 10.0;
   }

//====undo=====
   undo_vec.push_back(8);
   temp_data.push_back(undo_temp);
}

void ImgProc::rms_contrast()
{
   if( img_data == nullptr ){ return; }
   vector<float> undo_temp;
   undo_temp.resize(Nsize);
   //try to improve my code to be compatible
   //float mean_r, mean_g, mean_b;
   float mean = 0.0;
   //float dev_r, dev_g, dev_b;
   float dev = 0.0;
   // mean_r = 0.0;
   // mean_g = 0.0;
   // mean_b = 0.0;
   // dev_r = 0.0;
   // dev_g = 0.0;
   // dev_b = 0.0;

   for(int j=0; j<Nc; j++)
   {
     mean = 0.0;
     dev = 0.0;
     //for undo,stroe the data into memory
     for(long i=j; i<Nsize; i=i+Nc)
     {
       undo_temp[i] = img_data[i];
     }


     for( long i=j;i<Nsize;i = i+Nc )
     {
        mean += img_data[i];
     }
     mean /= (float)single_Nsize;

     for( long i=j;i<Nsize;i = i+Nc )
     {
        dev += (img_data[i] - mean) * (img_data[i] - mean);
     }
     dev /= (float)single_Nsize;
     dev = sqrt(dev);

     for( long i=j;i<Nsize;i = i+Nc )
     {
         img_data[i] = (img_data[i] - mean) / dev;
     }

   }

//===============================================================================

//====undo=====
   undo_vec.push_back(9);
   temp_data.push_back(undo_temp);
}

//==============part for the undo system===========
void ImgProc::undo_step()
{
   if (undo_vec.size() == 0)
   {
     cout << "No more steps to undo\n";
     return;
   }
   int udo_type = 0;
   //udo_type = undo_vec.back();
   udo_type = undo_vec[undo_vec.size()-1];
   //std::cout << "type:"    <<  udo_type <<"\n";
   if (udo_type == 0)
   {
     cout << "Something is not working I guess\n";
     return;
   }

   else if (udo_type == 1)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = img_data[i] / 1.1;
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
     //std::cout << "type"    <<  udo_type <<"\n";
   }

   else if (udo_type == 2)
   {
     if( img_data == nullptr ){ return; }
  #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = img_data[i] * 1.1;
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 3)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = img_data[i] - 0.05;
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 4)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = img_data[i] + 0.05;
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 5)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = pow(img_data[i], (1.0/0.9));
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 6)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = pow(img_data[i], (1.0/1.1));
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 7)
   {
     if( img_data == nullptr ){ return; }
     int t = temp_data.size()-1;
     if(t<0){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = temp_data[t][i];
     }
     temp_data.pop_back();
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 8)
   {
     if( img_data == nullptr ){ return; }
     int t = temp_data.size()-1;
     if(t<0){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = temp_data[t][i];
     }
     temp_data.pop_back();
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 9)
   {
     if( img_data == nullptr ){ return; }
     int t = temp_data.size()-1;
     if(t<0){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = temp_data[t][i];
     }
     temp_data.pop_back();
     undo_vec.pop_back();
     std::cout << "undo\n";
   }

   else if (udo_type == 10)
   {
     if( img_data == nullptr ){ return; }
   #pragma omp parallel for
     for( long i=0;i<Nsize;i++ )
     {
         img_data[i] = 1.0 - img_data[i];
     }
     undo_vec.pop_back();
     std::cout << "undo\n";
   }
}

//===================overview of the struct(included in C system)==================
// struct tm
// {
//   int tm_sec; //0~59, 60 if leap second23:59:60
//   int tm_min; //0-59
//   int tm_hour; //0-23
//   int tm_mday; //1-31
//   int tm_mon; //0-11
//   int tm_year; //Year - 1900
//   int tm_wday; // 0-6, Sunday = 0
//   int tm_yday; //0-365, 1 Jan = 0
//   int tm_isdst; //>0 summer time on; <0 unusable;
//   char  *tm_zone; //zone name
// }
//========================output system==========================
void ImgProc::output_exr()
{
  time_t tt = time(NULL);
  struct tm * t = localtime(&tt);
  int ty = t->tm_year + 1900;
  int tm = t->tm_mon + 1;
  string tmp = to_string(ty) + to_string(tm)+ to_string(t->tm_mday) + to_string(t->tm_hour)
   + to_string(t->tm_min)+ "_" + to_string(t->tm_sec);
   if(undo_vec.size() == 0)
   {tmp+="original.exr";}
   else
   {
     int udo_type = undo_vec[undo_vec.size()-1];
     if (udo_type == 1||udo_type == 2)
     {
       tmp += "brightness.exr";
     }

     else if (udo_type == 3||udo_type == 4)
     {
       tmp+="bias.exr";
     }

     else if (udo_type == 5||udo_type == 6)
     {
       tmp+="gamma.exr";
     }

     else if (udo_type == 7)
     {
       tmp+="grayscale.exr";
     }

     else if (udo_type == 8)
     {
       tmp+="quantize.exr";
     }

     else if (udo_type == 9)
     {
       tmp+="rms_contrast.exr";
     }

     else if (udo_type == 10)
     {
       tmp+="compliment.exr";
     }
  }

  //name the file
  const char *filename = tmp.c_str();
  //ttmp = ty + tm + "foo.exr";
  //const char *filename = tostring(ty) + tostring(tm) + "foo.exr";
  std::cout << "file_name:"  <<  tmp <<"\n";
  const int xres = Nx, yres = Ny;
  const int channels = Nc;
  //unsigned char pixels[xres*yres*channels];

  //std::unique_ptr<ImageOutput> out = ImageOutput::create (filename);
  auto out = ImageOutput::create (filename);
  if(! out)
      return;
  //ImageSpec spec (xres, yres, channels, TypeDesc::UINT8);
  ImageSpec spec (xres, yres, channels, TypeDesc::FLOAT);
  out->open (filename, spec);
  //out->write_image (TypeDesc::UINT8, pixels);
  out->write_image (TypeDesc::FLOAT, img_data);
  out->close ();
}

void ImgProc::output_jpeg()
{
  time_t tt = time(NULL);
  struct tm * t = localtime(&tt);
  int ty = t->tm_year + 1900;
  int tm = t->tm_mon + 1;
  string tmp = to_string(ty) + to_string(tm)+ to_string(t->tm_mday) + to_string(t->tm_hour)
   + to_string(t->tm_min)+ "_" + to_string(t->tm_sec);
  if(undo_vec.size() == 0)
  {tmp+="original.jpg";}
  else
  {
    int udo_type = undo_vec[undo_vec.size()-1];
    if (udo_type == 1||udo_type == 2)
    {
      tmp += "brightness.jpg";
    }

    else if (udo_type == 3||udo_type == 4)
    {
      tmp+="bias.jpg";
    }

    else if (udo_type == 5||udo_type == 6)
    {
      tmp+="gamma.jpg";
    }

    else if (udo_type == 7)
    {
      tmp+="grayscale.jpg";
    }

    else if (udo_type == 8)
    {
      tmp+="quantize.jpg";
    }

    else if (udo_type == 9)
    {
      tmp+="rms_contrast.jpg";
    }

    else if (udo_type == 10)
    {
      tmp+="compliment.jpg";
    }
  }


  //name the file
  const char *filename = tmp.c_str();
  std::cout << "file_name:"  <<  tmp <<"\n";
  const int xres = Nx, yres = Ny;
  const int channels = Nc;
  auto out = ImageOutput::create (filename);
  if(! out)
      return;
  ImageSpec spec (xres, yres, channels, TypeDesc::FLOAT);
  out->open (filename, spec);
  out->write_image (TypeDesc::FLOAT, img_data);
  out->close ();

}

//======================================================================================

long ImgProc::index(int i, int j, int c) const
{
   return (long) c + (long) Nc * index(i,j); // interleaved channels,rgbrgbrgb

   // return index(i,j) + (long)Nx * (long)Ny * (long)c; // sequential channels
}

long ImgProc::index(int i, int j) const
{
   return (long) i + (long)Nx * (long)j;
}



void img::swap(ImgProc& u, ImgProc& v)
{
   float* temp = v.img_data;
   int Nx = v.Nx;
   int Ny = v.Ny;
   int Nc = v.Nc;
   long Nsize = v.Nsize;

   v.Nx = u.Nx;
   v.Ny = u.Ny;
   v.Nc = u.Nc;
   v.Nsize = u.Nsize;
   v.img_data = u.img_data;

   u.Nx = Nx;
   u.Ny = Ny;
   u.Nc = Nc;
   u.Nsize = Nsize;
   u.img_data = temp;
}

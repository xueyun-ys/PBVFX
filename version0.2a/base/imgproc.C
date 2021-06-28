
#include <cmath>
#include "imgproc.h"
#include <string>
#include <iostream>

#include "Objects.h"
#include "operations.h"
#include "Camera.h"
#include "Color.h"
#include "Image.h"
#include "createFields.h"
#include "grids.h"

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

//===============================turntable=============================================
void ImgProc::turntable()
{
  Camera cam;
  float Pi = 3.1415926;
  float delta_cam = Pi/60.0;
  for(int i=112;i<120;i++)
  {

    cam.setEyeViewUp(Vector(25.0*sin(delta_cam*i), 0, 25.0*sin(Pi/2.0-delta_cam*i)), Vector(-25.0*sin(delta_cam*i), 0, -25.0*sin(Pi/2.0-delta_cam*i)), Vector(0, 1, 0));
    render(cam);
    output_jpeg();
    output_exr();
  }
}
//===============================turntable=============================================


//==================================================================================
//volume rendering
void ImgProc::render(Camera cam)//input: density field, color field
{
   if( img_data == nullptr ){ return; }
   //clear(Nx, Ny, 4);
   clear(1920, 1080, 4);
   //clear(360, 270, 4);
    //render avator
    Avator obj;
    ObjectDensityFieldAdd *dens_field = obj.getfield();
    ObjectColorFieldAdd *color_field = obj.getcolor();
    //obj.getfield()->eval()???????
    //std::cout<<"test4444"<<std::endl;
    //Grids gr(Vector(-10.0, -10.0, -10.0), 201,201,201,0.1);
    Grids *gr = new Grids(Vector(-6.0, -6.0, 6.0), 61,61,61,0.2);
    //std::cout<<"test5555"<<std::endl;
    FieldsToGrids *ftg = new FieldsToGrids(*gr, dens_field);
    //std::cout<<"test6666"<<std::endl;
    GridsToFields *gtf = new GridsToFields(ftg->get_grid());
    //std::cout<<"test7777"<<std::endl;
  	// using namespace std;
  	// using namespace lux;


    int test = 0;

  	//Camera cam;
  	//Image img;
    //#pragma omp parallel for  //can't use two for two loops like this
  	//outloop should be Ny!
    for (int v = 0; v< Ny ;v++)
  	{
      #pragma omp parallel for
  		for (int u = 0;u< Nx;u++)
  		{
  			Vector Np;
  			float T = 1.0;
  			float K = 0.9;
  			double del_s = 0.1;//0.01
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
          //std::cout<<"test3333"<<std::endl;
//1.3 5.3 1.8
  				//del_t = exp(-1.0 * K * del_s * gtf->eval(X));
          //Vector temp(1.3, 5.3, 1.8);
          //cout<<"ds_field_X:  "<<dens_field->eval(temp)<<endl;cout<<"gtf:  "<< gtf->eval(temp)<<endl;
          //do{cout<<'\n'<<"press any key to continue";}while(cin.get()!='\n');
          float dsf_eval = dens_field->eval(X);
          //if(dsf_eval!=0){cout<<"dsf_eval_X:  "<<X.X()<<" "<<X.Y()<<" "<<X.Z()<<endl;cout<<"gtf:  "<< gtf->eval(X)<<endl;}
          del_t = exp(-1.0 * K * del_s * dsf_eval);

          // if(dens_field->eval(X)!=0)
          // {
          //   std::cout<<"x: "<<X.X()<<"      y: "<<X.X()<<"      z: "<<X.X()<<std::endl;
          //
          // }
          //cout<<"1111111111111111111:    "<<endl;
          //double test = s2->eval(X).Z();
  				Lp = Lp + color_field->eval(X) * T * (1.0 - del_t) / K;
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
  				//Lp.set(Lp.X(), Lp.Y(), Lp.Z(), alpha);
          //if(Lp.Z()>0)
            //cout<<"!!!:    "<<Lp.Z()<<endl;
  			}
        long dex = index(u,v,0);
        //img_data[index(i,j,c)] = pixel[c];
        img_data[dex] = Lp.X();
        img_data[dex+1] = Lp.Y();
        img_data[dex+2] = Lp.Z();
        img_data[dex+3] = alpha;//Lp.W();
        //test ++;
        //if(Lp.Z()>0)
          //cout<<"!!!:    "<<Lp.Z()<<endl;
        //cout<<"test:    "<<test<<endl;
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

void ImgProc::render()
{
   if( img_data == nullptr ){ return; }
    //render avator

    //ConstantField *densityVolume = new ConstantField(1.0);
  	ObjectSphere *densityVolume = new ObjectSphere(10);
    //ObjectBox *densityVolume = new ObjectBox(1.0, 5.0);
    Color c(0.5d,0.0d,0.5d,0.5d);
  	ConstantColor *colorVolume = new ConstantColor(c);
    //colorVolume->set(0.0, 0.0, 0.5, 0.0);
    ObjectMask *s1 = new ObjectMask(densityVolume);
    ObjectColorDensityFieldMultiple *s2 = new ObjectColorDensityFieldMultiple(s1, colorVolume);
  	//create objects

    int test = 0;

  	Camera cam;
  	//Image img;
    //#pragma omp parallel for
  	//outloop should be Ny!
    for (int v = 0; v< Ny ;v++)
  	{
      //#pragma omp parallel for
  		for (int u = 0;u< Nx;u++)
  		{
  			Vector Np;
  			float T = 1.0;
  			float K = 0.9;
  			double del_s = 0.5;
  			double s;
  			float alpha;
  			float del_t;
        double x = (double)u / (double)Nx;
        double y = (double)v / (double)Ny;
  			Np = cam.view(x, y);
  			Vector X = cam.eye() + Np * cam.nearPlane();
  			s = cam.nearPlane();
  			 //Lp = (0, 0, 0, 0)
  			Color Lp;


  			while (s < cam.farPlane() && T > 0.000001)
  			{
  				X = X + del_s * Np;
  				s += del_s;
  				del_t = exp(-1.0 * K * del_s * densityVolume->eval(X));
  				Lp = Lp + s2->eval(X) * T * (1.0 - del_t) / K;
  				T *= del_t;
          alpha = 1.0 - T;
  				Lp.set(Lp.X(), Lp.Y(), Lp.Z(), alpha);
  			}
        long dex = index(u,v,0);
        img_data[dex] = Lp.X();
        img_data[dex+1] = Lp.Y();
        img_data[dex+2] = Lp.Z();
        //img_data[dex+3] = alpha;
        test ++;

  		}

    }
}















// #include "Objects.h"
// #include "Camera.h"
// #include "Color.h"
// #include "Image.h"
//
// int main(int argc, char* argv[])
// {
// 	using namespace std;
// 	using namespace lux;
//
// 	Volume<float> densityVolume;
// 	Volume<Color> colorVolume;
// 	//create objects
//
// 	Camera cam;
// 	Image img;
// 	for (int u = 0;u< ???;u++)
// 	{
// 		for (int v = 0; v<??? ;v++)
// 		{
// 			Vector Np;
// 			float T = 1.0;
// 			float K = 0.9;
// 			double del_s = 0.05;
// 			double s;
// 			float alpha;
// 			float del_t;
// 			//cam.setEyeViewUp(Vector(0, -0, -2.0), Vector(0, 0, 1.0), Vector(0, 0, -1));
// 			//cam.setRight();
// 			//cam.setAspectRatio(10);
// 			//cam.setFov(178);
// 			Np = cam.view(x, y);
// 			//Vector Xp = cam.eye() + u * cam.right() + v * cam.up();
// 			//Vector Np = (Xp - cam.eye()) / (cam.view() * (Xp - cam.eye())) - cam.view();
// 			Vector X = cam.eye() + Np * cam.nearPlane();
// 			s = cam.nearPlane();
// 			 //Lp = (0, 0, 0, 0)
// 			Color Lp;
//
// 			while (s < cam.farPlane() && T > 0.000001)
// 			{
// 				X = X + del_s * Np;
// 				s += del_s;
// 				del_t = exp(-1.0 * K * del_s * densityVolume.eval(X));
// 				Lp = Lp + colorVolume.eval(X) * T * (1 - del_t) / K;
// 				alpha = 1.0 - T;
// 				T *= del_t;
// 				Lp.set(Lp.X(), Lp.Y(), Lp.Z(), alpha);
// 			}
// 		}
// 	}
//
//
// 	return 0;
// }

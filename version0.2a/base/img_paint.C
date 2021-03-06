//------------------------------------------------
//
//  img_paint
//
//
//-------------------------------------------------




#include <cmath>
#include <omp.h>
#include "imgproc.h"
#include "CmdLineFind.h"
#include <vector>
#include "Camera.h"


#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.


#include <iostream>
#include <stack>


using namespace std;
using namespace img;

ImgProc image;




void setNbCores( int nb )
{
   omp_set_num_threads( nb );
}

void cbMotion( int x, int y )
{
}

void cbMouse( int button, int state, int x, int y )
{
}

void cbDisplay( void )
{
   glClear(GL_COLOR_BUFFER_BIT );
   if(image.depth() <= 3)
      glDrawPixels( image.nx(), image.ny(), GL_RGB, GL_FLOAT, image.raw() );
   else
      glDrawPixels( image.nx(), image.ny(), GL_RGBA, GL_FLOAT, image.raw() );
   glutSwapBuffers();
}

void cbIdle()
{
   glutPostRedisplay();
}

void cbOnKeyboard( unsigned char key, int x, int y )
{
   switch (key)
   {
      case 'r':
   image.turntable();
   cout << "render avator\n";
   break;

   //    case 'R':
   // image.render();
   // cout << "render avator\n";
   // break;

      case 'c':
	 image.compliment();
	 cout << "Compliment\n";
	 break;

      case 'C':
   image.rms_contrast();
   cout << "rms_Contrast\n";
   break;

      case 'V':
	 image.brightness_up();
	 cout << "Increase Brightness\n";
	 break;

      case 'v':
	 image.brightness_down();
	 cout << "Decrease Brightness\n";
	 break;

      case 'B':
   image.bias_up();
   cout << "Increase Bias\n";
   break;

      case 'b':
	 image.bias_down();
	 cout << "Decrease Bias\n";
	 break;

      case 'G':
	 image.gamma_up();
	 cout << "Increase Gamma\n";
	 break;

      case 'g':
	 image.gamma_down();
	 cout << "Decrease Gamma\n";
	 break;

      case 'w':
	 image.grayscale();
	 cout << "Create Grayscale\n";
	 break;

      case 'q':
   image.quantize();
   cout << "Quantize\n";
   break;

      case 'u':
   image.undo_step();
   //cout << "undo\n";
   break;

      case 'o':
   image.output_exr();
   cout << "output the exr image into memory\n";
   break;

      case 'j':
   image.output_jpeg();
   cout << "output the jpeg image into memory\n";
   break;
   }
}

void PrintUsage()
{
   cout << "img_paint keyboard choices\n";
   cout << "r         render avator\n";
   cout << "c         compliment\n";
   cout << "C         rms_contrast\n";
   cout << "v         brightness_down\n";
   cout << "V         brightness_up\n";
   cout << "b         bias_down\n";
   cout << "B         bias_up\n";
   cout << "g         gamma_down\n";
   cout << "G         gamma_up\n";
   cout << "w         grayscale\n";
   cout << "q         quantize\n";
   cout << "u         undo_step\n";
   cout << "o         output exr file\n";
   cout << "j         output jpg file\n";
}


int main(int argc, char** argv)
{
   lux::CmdLineFind clf( argc, argv );

   setNbCores(8);

   string imagename = clf.find("-image", "", "Image to drive color");

   clf.usage("-h");
   clf.printFinds();
   PrintUsage();

   image.load(imagename);


   // GLUT routines
   glutInit(&argc, argv);

   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   glutInitWindowSize( image.nx(), image.ny() );

   // Open a window
   char title[] = "img_paint";
   glutCreateWindow( title );

   glClearColor( 1,1,1,1 );

   glutDisplayFunc(&cbDisplay);
   glutIdleFunc(&cbIdle);
   glutKeyboardFunc(&cbOnKeyboard);
   glutMouseFunc( &cbMouse );
   glutMotionFunc( &cbMotion );

   glutMainLoop();
   return 1;
};

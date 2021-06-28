//*******************************************************************
//
//   CreateOBJ.h
//
//
//
//
//
//*******************************************************************

#ifndef __CREATEOBJ_H__
#define __CREATEOBJ_H__

#include <cmath>
#include <string>
#include <iostream>

//#include "Objects.h"
//#include "operations.h"
//#include "Color.h"
//#include "Image.h"
#include <fstream>
#include <stdlib.h>
#include <sstream>
//stdlib.h里面定义了五种类型、一些宏和通用工具函数。类型例如size_t、wchar_t、div_t、ldiv_t和lldiv_t；宏例如EXIT_FAILURE、EXIT_SUCCESS、RAND_MAX和MB_CUR_MAX等等；常用的函数如malloc()、calloc()、realloc()、free()、system()、atoi()、atol()、rand()、 srand()、exit()等等。
#include <string>
#include <iostream>


//using namespace img;
using namespace std;
//using namespace lux;

class OBJ_model
{
  public:
    OBJ_model(string path)
    {

    }
    OBJ_model()
    {
      //======================================
      	string s, str, s1, s2, s3, s4;
      	ifstream inf;
      	inf.open("bunny.obj");
      	int vn = 0;
      	int vnum = 0;
      	int fnum = 0;


      	while (getline(inf, s))
      	{
      		istringstream in(s);
      		in >> s1 >> s2 >> s3 >> s4;
      		if (s[0]=='v')
      		{
      			if (s[1] == 'n')
      			{
      				vn++;
              cout<<"wrong!"<<s[0]<<"  "<<s[1]<<"  "<<s[2]<<endl;
      			}
      			else
      			{
      				vnum++;
      			}

      		}
    //-------------------------------------------------------
          if(s[0]=='v'){
              if(s[1]=='t')
              {
              }
              else if(s[1]=='n')
              {
              }
              else{
                  istringstream in(s);
                  vector<double> p;
                  double x, y, z;
                  double scale_up = 100;
                  string head;
                  in>>head>>x>>y>>z;
                  //cout <<"x: "<< x << "y: "<< y <<"z: "<< z <<endl;
                  //do{cout<<'\n'<<"press any key to continue!!!";}while(cin.get()!='\n');
                  p.push_back(x*scale_up);
                  p.push_back(y*scale_up-10.0);
                  p.push_back(z*scale_up);
                  vp.push_back(p);
              }
          }

    //-------------------------------------------------------------------------------
      		float a = 0, b = 0, c = 0;
      		if (s[0] == 'f')
      		{
            vector<int> fi;
      			for (int k = 0; k<s2.size() ; k++)//s2[k] != ' '
      			{
      				a = a * 10 + (s2[k] - 48);
      			}
            fi.push_back(a);
    //do{cout<<'\n'<<"press any key to continue!!!";}while(cin.get()!='\n');
      			for (int k = 0; k<s3.size(); k++)
      			{
      				b = b * 10 + (s3[k] - 48);
      			}
            fi.push_back(b);
      			for (int k = 0; k<s4.size(); k++)
      			{
      				c = c * 10 + (s4[k] - 48);
      			}
            fi.push_back(c);
            // cout <<"a: "<< a << "b: "<< b <<"c: "<< c <<endl;
            // do{cout<<'\n'<<"press any key to continue!!!";}while(cin.get()!='\n');
            f_index.push_back(fi);
      			fnum++;

      		}
      	}

      	inf.close();

    //============================================

   }

   ~OBJ_model(){}

   vector<vector<int>> getface()
   {
      return f_index;
   }
   vector<vector<double>> getpoints()
   {
      return vp;
   }

  private:
    vector<vector<double>> vp;
    vector<vector<int>> f_index;
};



#endif

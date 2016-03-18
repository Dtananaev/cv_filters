/*
 * nlm.cpp
 *
 *  Created on: February 12, 2015
 *      Author: Denis Tananaev
 */


#include "CMatrix.h"
#include "CFilter.h"
#include "NMath.h"
#include <iostream>
#include <cstdio>
#include <ctime>
#include <omp.h>//!!!!uncomment in case if you have installed omp libraries!!!
#include "CTensor.h"


inline float patch_distance( const CMatrix<float>& img, 
                             int x1, int y1, 
                             int x2, int y2, 
                             int patch_radius,
                             float* gauss_lut_center )
{
  const int x_size = img.xSize();
  const int y_size = img.ySize();

  float ssd = 0;
  for( int ty = -patch_radius; ty <= patch_radius; ++ty )
  for( int tx = -patch_radius; tx <= patch_radius; ++tx )
  {
    // clamp coordinates
    int p1x = std::min(x_size-1,std::max(0,x1+tx));//start and end of patch
    int p1y = std::min(y_size-1,std::max(0,y1+ty));
    int p2x = std::min(x_size-1,std::max(0,x2+tx));
    int p2y = std::min(y_size-1,std::max(0,y2+ty));
    float tmp = img(p1x,p1y)-img(p2x,p2y);
    float gauss_w = *(gauss_lut_center+tx);

    gauss_w *= *(gauss_lut_center+ty);
		  
			
  ssd += tmp*tmp*gauss_w;
	
  }


  return ssd;
}

CMatrix<float> nlmfilter(CMatrix<float> aImage, int window_radius, int patch_radius,float sigma){
	//initialize constant
const float inv_sqr_sigma = 1/(sigma*sigma);//1/sigma^2
  const int x_size = aImage.xSize();
  const int y_size = aImage.ySize();
  CMatrix<float> result(x_size,y_size);
   
  // create a gauss lut for the function patch_distance()
  float* gauss_lut = new float[2*patch_radius+1];
  float* gauss_lut_center = gauss_lut+patch_radius;
  for( int x = -patch_radius; x <= patch_radius; ++x )
    *(gauss_lut_center+x) = exp(-0.5*x*x/(patch_radius*patch_radius));
 
#pragma omp parallel for//!!!!uncomment in case if you have installed omp libraries!!!
   for( int y = 0; y < y_size; ++y )
  for( int x = 0; x < x_size; ++x )
	{
	const int x1 = std::max(0,x-window_radius);
    const int y1 = std::max(0,y-window_radius);
    const int x2 = std::min(x_size-1,x+window_radius);
    const int y2 = std::min(y_size-1,y+window_radius);
	

   float sum = 0;
    float new_value = 0;

    for( int ny = y1; ny <= y2; ++ny )
    for( int nx = x1; nx <= x2; ++nx )
    {
      float dsqr = patch_distance(aImage,x,y,nx,ny,patch_radius,gauss_lut_center);
      float w = exp(-dsqr*inv_sqr_sigma);
      new_value += w*aImage(nx,ny);
      sum += w;
    }
    result(x,y) = new_value/sum;
  }

  return result;


}




int main(int argc, char** argv){
	

 std::string fileNameInput;
  int window_radius=10;
  int patch_radius=3;
    float sigma=100;

    
    if (argc==2){
        fileNameInput=argv[1];

    } else if (argc==3){
       fileNameInput=argv[1];
        patch_radius=atoi(argv[2]);
 
    } else if (argc==4){
         fileNameInput=argv[1];
        patch_radius=atoi(argv[2]);
        window_radius=atoi(argv[3]);
     } else if (argc==5){
         fileNameInput=argv[1];
        patch_radius=atoi(argv[2]);
        window_radius=atoi(argv[3]);  
        sigma=atoi(argv[4]); 
    }else{
        std::cout<<"!!!WRONG INPUT!!!"<<"\n";
        std::cout<<"Usage: nlm inputfile  <path radius> <window radius> <sigma>"<<"\n";
        std::cout<<"The command should contain at least input file name. The default nlm path_radius=3, window=10, sigma=100."<<"\n";
        return 0;    
    }

	CTensor<float> image_color;
	image_color.readFromPPM((fileNameInput+".ppm").c_str());
	
	CMatrix<float> red(image_color.xSize(),image_color.ySize()), filter_red(image_color.xSize(),image_color.ySize());
    CMatrix<float> green(image_color.xSize(),image_color.ySize()), filter_green(image_color.xSize(),image_color.ySize());
    CMatrix<float> blue(image_color.xSize(),image_color.ySize()), filter_blue(image_color.xSize(),image_color.ySize());
  
    CTensor<float> result_color(image_color.xSize(),image_color.ySize(),image_color.zSize(),0);
    CMatrix<float> Result(image_color.xSize(),image_color.ySize());
    std::cout<<"filter works wait... it could took 2-5 minutes "<<"\n";
	    clock_t start;
        double duration;

    start = std::clock();
	filter_red=nlmfilter(image_color.getMatrix(0),window_radius,patch_radius,sigma);
     std::cout<<"1/3 Done. "<<"\n";
	filter_green=nlmfilter(image_color.getMatrix(1),window_radius,patch_radius,sigma);
         std::cout<<"2/3 Done. "<<"\n";
	filter_blue=nlmfilter(image_color.getMatrix(2),window_radius,patch_radius,sigma);
             std::cout<<"3/3 Done. "<<"\n";
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

   std::cout<<"Time of nlm filter: "<< duration<<" seconds"<<'\n';
   
    result_color.putMatrix(filter_red,0);
    result_color.putMatrix(filter_green,1);
    result_color.putMatrix(filter_blue,2);
	result_color.writeToPPM((fileNameInput+"_nlm.ppm").c_str());

	return 0;
}

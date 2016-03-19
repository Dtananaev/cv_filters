/*
 * wavelet.cpp
 *
 *  Created on: November 25, 2014
 *      Author: Denis Tananaev
 */

#include "CMatrix.h"
#include <stdlib.h> 
#include <ctime>
#include <cmath>
#include <string>
#include "CTensor.h"
#include <string>




//Back transform
CMatrix<float> backtransform( CMatrix<float> aImage,int level)	
{
	int n=pow(2.0,level);

	int Xsize=aImage.xSize()/n;
 
	int Ysize=aImage.ySize()/n;

	CMatrix<float> CK(Xsize, Ysize,255);
	CMatrix<float> DV(Xsize, Ysize,255);
	CMatrix<float> DD(Xsize, Ysize,255);
	CMatrix<float> DH(Xsize, Ysize,255);
	CMatrix<float> C(2*Xsize,2*Ysize,255);
	
	for (int y = 0; y < CK.ySize(); y++)
	for (int x = 0; x < CK.xSize(); x++){
	 CK(x,y)= aImage(x,y);
	 DH(x,y)=aImage(x+CK.xSize(),y);
	DV(x,y)= aImage(x,y + CK.ySize());
	DD(x,y)=aImage(x+CK.xSize(),y + CK.ySize());
	}

	 for (int y = 0; y < CK.ySize(); y++)
    for (int x = 0; x < CK.xSize(); x++) 
	{
		C(2*x,2*y)=CK(x,y)+DH(x,y)+DV(x,y)+DD(x,y);
		C(2*x+1,2*y)=CK(x,y)-DH(x,y)+DV(x,y)-DD(x,y);
		C(2*x,2*y+1)=CK(x,y)+DH(x,y)-DV(x,y)-DD(x,y);
		C(2*x+1,2*y+1)=CK(x,y)-DH(x,y)-DV(x,y)+DD(x,y);
    }
	 for (int y = 0; y < C.ySize(); y++)
    for (int x = 0; x < C.xSize(); x++) 
	{
	aImage(x,y)=C(x,y);
	}

	if (level>1)
	{aImage=backtransform( aImage, --level);}

	return aImage;
}

// Compute wavelet coefficients multi-scale decomposition
CMatrix<float> wavetransform( CMatrix<float> aImage, int level)	
{
	int Xsize=aImage.xSize();
	int Ysize=aImage.ySize();
	CMatrix<float> C(aImage);
	CMatrix<float> CK(Xsize/2, Ysize/2,255);
	CMatrix<float> DV(Xsize/2, Ysize/2,255);
	CMatrix<float> DD(Xsize/2, Ysize/2,255);
	CMatrix<float> DH(Xsize/2, Ysize/2,255);
	CMatrix<float> Result(aImage.xSize(),aImage.ySize());
 for (int y = 0; y < CK.ySize(); y++)
    for (int x = 0; x < CK.xSize(); x++) {
		CK(x,y)=0.25*(C(2*x,2*y)+C(2*x+1,2*y)+C(2*x,2*y+1)+C(2*x+1,2*y+1));

		DH(x,y)=0.25*(C(2*x,2*y)+C(2*x,2*y+1)-C(2*x+1,2*y)-C(2*x+1,2*y+1));

		DV(x,y)=0.25*(C(2*x,2*y)+C(2*x+1,2*y)-C(2*x,2*y+1)-C(2*x+1,2*y+1));

		DD(x,y)=0.25*(C(2*x,2*y)-C(2*x+1,2*y)-C(2*x,2*y+1)+C(2*x+1,2*y+1));

    }
	C=CK;
	if (level>1)
		CK=wavetransform( CK, --level);

	for (int y = 0; y < CK.ySize(); y++)
	for (int x = 0; x < CK.xSize(); x++){
	  Result(x,y) = CK(x,y);
	  Result(x+CK.xSize(),y) = DH(x,y);
	  Result(x,y + CK.ySize()) = DV(x,y);
	  Result(x+CK.xSize(),y + CK.ySize()) = DD(x,y);
	}

	return Result;

}


///For colored images
void image2rgb(CTensor<float> image_color, CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue){
    for(int x=0;x<image_color.xSize();x++)
	 for(int y=0;y<image_color.ySize();y++)
	 {
	 red(x,y)=image_color(x,y,0);
	 green(x,y)=image_color(x,y,1);
	 blue(x,y)=image_color(x,y,2);
	 }

}

void rgb2image(CTensor<float>& result_color, CMatrix<float> red, CMatrix<float> green, CMatrix<float> blue){
	for(int x=0;x<result_color.xSize();x++)
	for(int y=0;y<result_color.ySize();y++)
	 {
	result_color(x,y,0)=red(x,y);
	result_color(x,y,1)=green(x,y);
	result_color(x,y,2)=blue(x,y);
	 }
}

void applyWaveDecompose(CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, int level){

     red=wavetransform(red, level);	
     green=wavetransform(green, level);	
    blue=wavetransform(blue, level);	

}

void applyWaveRecover(CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, int level){

     red=backtransform(red, level);	
     green=backtransform(green, level);	
    blue=backtransform(blue, level);	

}


int main(int argc, char** argv) {

    std::string fileNameInput;
    int level=1;
    
    if (argc==2){
        fileNameInput=argv[1];

    } else if (argc==3){
       fileNameInput=argv[1];
        level=atoi(argv[2]);
 
 
    }else{
        std::cout<<"!!!WRONG INPUT!!!"<<"\n";
        std::cout<<"Usage: wavelet inputfile  <level of decomposition> "<<"\n";
        std::cout<<"The command should contain at least input file name. The default level of wavelet decomposition 1."<<"\n";
        return 0;    
    }

    std::cout<<"file name is "<<fileNameInput.c_str()<<"\n";

    //Color filter

    CTensor<float> image_color;

    image_color.readFromPPM((fileNameInput+".ppm").c_str());

    CMatrix<float> red(image_color.xSize(),image_color.ySize()),filter_red(image_color.xSize(),image_color.ySize());
    CMatrix<float> green(image_color.xSize(),image_color.ySize()),filter_green(image_color.xSize(),image_color.ySize());
    CMatrix<float> blue(image_color.xSize(),image_color.ySize()),filter_blue(image_color.xSize(),image_color.ySize());

    CTensor<float> result_color(image_color.xSize(),image_color.ySize(),image_color.zSize(),0);

    image2rgb(image_color,red, green, blue);

    applyWaveDecompose(red, green, blue,level);

    rgb2image(result_color,red, green, blue);

    result_color.writeToPPM((fileNameInput+"_coeff.ppm").c_str());

    applyWaveRecover(red, green, blue,level);
    
    rgb2image(result_color,red, green, blue);
  
   result_color.writeToPPM((fileNameInput+"_recovered.ppm").c_str());

  return 0;
}

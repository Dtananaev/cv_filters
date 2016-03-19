/*
* gaussn.cpp
*
*  Created on: November 11, 2014
*      Author: Denis Tananaev
*/

#include "CMatrix.h"
#include <stdlib.h> 
#include <ctime>
#include <cmath>
#include <string>
#include "CTensor.h"
#include <string>




//add Gaussian noise
CMatrix<float> addnoiseg(CMatrix<float> aImage, float mu, float var)
{
	srand(time(NULL));
	CMatrix<float> C(aImage);
	double U,V, N, M, n, m;	
	// add noise
	for (int y = 0; y < C.ySize()-1; y++)
		for (int x = 0; x < C.xSize()-1; x++) {
			double U=rand()/((double)RAND_MAX);
			double V=rand()/((double)RAND_MAX);

			//U can't be 0 because of log
			if(U==0){U=U+0.1;}
			// N and M normal random values N(0, 1)
			N=sqrt(-2.0*log(U)) * cos(2.0*3.1415926536*V);
			M=sqrt(-2.0*log(U)) * sin(2.0*3.1415926536*V);
			// n and m normap random variables N(mean, variance)
			n=mu+sqrt(var*var)*N;
			m=mu+sqrt(var*var)*M;

			aImage(x,y)=aImage(x,y)+n;
			if (aImage(x,y) > 255){
				aImage(x,y) = 255;}
			else if (aImage(x,y)< 0){
				aImage(x,y) = 0;}
			else
			{	 aImage(x,y)= aImage(x,y);}
		}
		return aImage;	
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

void addNoise( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, float mu, float var){
	red=addnoiseg(red,mu,var); 
	green=addnoiseg(green,mu,var); 
	blue=addnoiseg(blue,mu,var); 
}


int main(int argc, char** argv) {

	std::string fileNameInput;
	float mu=0;
	float sigma =30;

	if (argc==2){
		fileNameInput=argv[1];

	} else if (argc==3){
		fileNameInput=argv[1];
		mu=atoi(argv[2]);

	} else if (argc==4){
		fileNameInput=argv[1];
		mu=atoi(argv[2]);
		sigma=atoi(argv[3]);
	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: gaussn inputfile  <mean value> <variance>"<<"\n";
		std::cout<<"The command should contain at least input file name. The default noise parameters mean=0 and var=30."<<"\n";
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

	addNoise(red, green, blue,0,30);

	rgb2image(result_color,red, green, blue);

	result_color.writeToPPM((fileNameInput+"_gaussn.ppm").c_str());

	return 0;
}

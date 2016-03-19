/*
* boxf.cpp
*
*  Created on: November 23, 2014
*      Author: Denis Tananaev
*/

#include "CMatrix.h"
#include <stdlib.h> 
#include <ctime>
#include <cmath>
#include <string>
#include "CTensor.h"



//Cut Neumann boundaries
CMatrix<float> cut(CMatrix<float>& image,int border_size){ 

	CMatrix<float> realimage(image.xSize()-2*border_size,image.ySize()-2*border_size);
	for(int x=0;x<realimage.xSize();x++)
		for(int y=0; y<realimage.ySize();y++)
		{
			realimage(x,y)=image(x+border_size,y+border_size);
		}

		return realimage;
}


// Neumann boundry conditions
CMatrix<float> Neumann_bound(CMatrix<float> aImage,int border_size){

	CMatrix<float> result(aImage.xSize()+border_size*2,aImage.ySize()+border_size*2);
	result.fill(0);
	//center matrix
	for(int x=0;x<aImage.xSize();x++)
		for(int y=0;y<aImage.ySize();y++)
		{
			result(x+border_size,y+border_size)=aImage(x,y);
		}
		//Top
		for(int x=0;x<aImage.xSize();x++)
			for(int y=0;y<border_size;y++)
			{
				result(x+border_size,y)=aImage(x,border_size-1-y);
			}
			//Bottom
			for(int x=0;x<aImage.xSize();x++)
				for(int y=0;y<border_size;y++)
				{
					result(x+border_size,y+aImage.ySize()+border_size)=aImage(x,aImage.ySize()-1-y);
				}
				//left side
				for(int x=0;x<border_size;x++)
					for(int y=0;y<aImage.ySize();y++)
					{
						result(x,y+border_size)=aImage(border_size-1-x,y);
					}

					//right side
					for(int x=0;x<border_size;x++)
						for(int y=0;y<aImage.ySize();y++)
						{
							result(x+aImage.xSize()+border_size,y+border_size)=aImage(aImage.xSize()-1-x,y);
						}
						//up left square
						for(int x=0;x<border_size;x++)
							for(int y=0;y<border_size;y++)
							{
								result(x,y)=aImage(0,0);
							}
							//up right square
							for(int x=aImage.xSize()-1;x<(aImage.xSize()+border_size);x++)
								for(int y=0;y<border_size;y++)
								{
									result(x+border_size,y)=aImage(aImage.xSize()-1,0);
								}
								//down left square
								for(int x=0;x<border_size;x++)
									for(int y=aImage.ySize()-1;y<(aImage.ySize()+border_size);y++)
									{
										result(x,y+border_size)=aImage(0,aImage.ySize()-1);
									}
									//down right square
									for(int x=aImage.xSize()-1;x<(aImage.xSize()+border_size);x++)
										for(int y=aImage.ySize()-1;y<(aImage.ySize()+border_size);y++)
										{
											result(x+border_size,y+border_size)=aImage(aImage.xSize()-1,aImage.ySize()-1);
										}		
										return result;
}


//Box 2D
CMatrix<float> Box(float sigma)
{
	float sum=0;
	size_t filter_size = size_t(2*sigma+1);
	std::cout<<"filter_size_box= "<<" pixels"<<filter_size<<std::endl;
	CMatrix<float> filter(filter_size,filter_size);
	if ( filter_size % 2 == 0 )
	{ ++filter_size;}
	filter.fill(1/(filter_size*(2*sigma+1)));
	return filter;

}

//Box filter
CMatrix<float> Bfilter(CMatrix<float> Box, CMatrix<float> boundary_Image,float sigma,int iteration)
{
	CMatrix<float> image=boundary_Image;
	//CMatrix<float> result( boundary_Image.xSize(),boundary_Image.ySize());
	int border=sigma;
	for(int x=border;x<boundary_Image.xSize()-border;x++)
		for(int y=border;y<boundary_Image.ySize()-border;y++)
		{
			image(x,y)=0;
		}
		//Filtering

		for(int x=border;x<boundary_Image.xSize()-border;x++)
			for(int y=border;y<boundary_Image.ySize()-border;y++)
			{
				for(int i=0;i<Box.xSize();i++)
					for (int j=0;j<Box.ySize();j++)
					{
						image(x,y)+=boundary_Image(x-border+i,y-border+j)*Box(i,j);
					}


			}
			boundary_Image=image;
			if(iteration>1)
			{
				boundary_Image=Bfilter(Box,boundary_Image,sigma,--iteration);
			}
			return boundary_Image;

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


void applyBoxf( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, int  sigma, int iteration){

	//apply Neumann boundary conditions
	int border_size=sigma;

	red=Neumann_bound(red,border_size);
	green=Neumann_bound(green,border_size);
	blue=Neumann_bound(blue,border_size);
	//Create Gaussian kernel and apply Gauss filter
	CMatrix<float> kernel=Box(sigma);

	red=Bfilter(kernel, red, sigma, iteration);
	green=Bfilter(kernel,green,sigma, iteration);
	blue=Bfilter(kernel, blue, sigma, iteration);

	//cut the Neumann boundaries
	red=cut(red, border_size);
	green=cut(green, border_size);
	blue=cut(blue, border_size);
}


int main(int argc, char** argv) {

	std::string fileNameInput;
	int sigma=3;
	int iteration=1; 


	if (argc==2){
		fileNameInput=argv[1];

	} else if (argc==3){
		fileNameInput=argv[1];
		sigma=atoi(argv[2]);
	} else if (argc==4){
		fileNameInput=argv[1];
		sigma=atoi(argv[2]);
		iteration=atoi(argv[3]);
	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: boxf inputfile <size of box> <number of iterations of  box filter>"<<"\n";
		std::cout<<"The command should contain at least input file name. The default  kernel size sigma=3 and iteration=1."<<"\n";
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

	applyBoxf(red, green, blue, sigma, iteration);

	rgb2image(result_color,red, green, blue);

	result_color.writeToPPM((fileNameInput+"_boxf.ppm").c_str());



	return 0;
}

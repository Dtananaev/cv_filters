/*
* diffusion.cpp
*
*  Created on: December 24, 2014
*      Author: Denis Tananaev
*/

#include "CFilter.h"
#include "CMatrix.h"
#include "CVector.h"
#include "CTensor.h"
#include "CTensor4D.h"
#include "NMath.h"
#include <iostream>

using namespace std;



//Cut Neumann boundaries
void cut(CMatrix<float>& image,int border_size)
{ 
	CMatrix<float> realimage(image.xSize()-2*border_size,image.ySize()-2*border_size);
	for(int x=0;x<realimage.xSize();x++)
		for(int y=0; y<realimage.ySize();y++)
		{
			realimage(x,y)=image(x+border_size,y+border_size);
		}
		image=image(image.xSize()-2*border_size,image.ySize()-2*border_size);
		image=realimage;
		//return realimage;
}
// Neumann boundry conditions
CMatrix<float> Neumann_bound(CMatrix<float> aImage,int border_size)
{
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

//PSNR function
double PSNR(CMatrix<float> aImage, CMatrix<float> noisyImage)
{
	double psnr,dif=0;	
	for(int i=0;i<aImage.xSize();i++)
		for(int j=0;j<aImage.ySize();j++)
		{dif+=pow(aImage(i,j)-noisyImage(i,j),2);}
		psnr=10*log10((pow(aImage.max()-aImage.min(),2)*aImage.size())/dif);
		return psnr;
}

// nonlinear diffusion dual algorithm
CMatrix<float> DUALnlindiffusion(CMatrix<float> aImage, float t,float alpha, int border, int k){

	CMatrix<float> result(aImage);
	CMatrix<float> p(aImage.xSize(),aImage.ySize(),0);

	CMatrix<float> dx_p(p); // p1
	CMatrix<float> dy_p(p); //p2
	CMatrix<float> p1(p); // ph1
	CMatrix<float> p2(p); // ph2

	CMatrix<float> z(p);

	for(int i=0;i<k;++i)
	{
		for(int x=border;x<aImage.xSize()-border;x++)
			for(int y=border;y<aImage.ySize()-border;y++)
			{
				{
					//backward difference
					z(x,y)=dx_p(x,y)-dx_p(x-1,y)+dy_p(x,y)-dy_p(x,y-1);
					//forward difference


					p1(x,y)=dx_p(x,y)+(t/alpha)*((aImage(x+1,y)+alpha*z(x+1,y))-(aImage(x,y)+alpha*z(x,y)));
					p2(x,y)=dy_p(x,y)+(t/alpha)*((aImage(x,y+1)+alpha*z(x,y+1))-(aImage(x,y)+alpha*z(x,y)));


					float norm=sqrtf(p1(x,y)*p1(x,y)+p2(x,y)*p2(x,y));
					if (norm>1)
					{dx_p(x,y)=p1(x,y)/norm;
					dy_p(x,y)=p2(x,y)/norm;
					}
					else
					{
						dx_p(x,y)=p1(x,y);
						dy_p(x,y)=p2(x,y);
					}
				}

			}
	}

	for(int x=border;x<aImage.xSize()-border;x++)
		for(int y=border;y<aImage.ySize()-border;y++)
		{
			result(x,y)=aImage(x,y)+alpha*(dx_p(x,y)-dx_p(x-1,y)+dy_p(x,y)-dy_p(x,y-1));
		}
		return result;
}

// nonlinear diffusion TV flow algorithm

CMatrix<float> TVnlindiffusion(CMatrix<float> aImage, float t,float e, int border, int k){

	CDerivative<float> aDerivative(3);
	CMatrix<float> dx(aImage.xSize(),aImage.ySize());
	CMatrix<float> dy(aImage.xSize(),aImage.ySize());
	NFilter::filter(aImage,dx,aDerivative,1);
	NFilter::filter(aImage,dy,1,aDerivative);
	CMatrix<float> result(aImage.xSize(),aImage.ySize(),255);


	for(int x=border;x<aImage.xSize()-border;x++)
		for(int y=border;y<aImage.ySize()-border;y++)
		{

			double g_1=0.5*(1/sqrt(pow(dx(x+1,y),2)+pow(dy(x+1,y),2)+pow(e,2))+ 1/sqrt(pow(dx(x,y),2)+pow(dy(x,y),2)+pow(e,2)));

			double g_2=0.5*(1/sqrt(pow(dx(x-1,y),2)+pow(dy(x-1,y),2)+pow(e,2))+ 1/sqrt(pow(dx(x,y),2)+pow(dy(x,y),2)+pow(e,2)));

			double g_3=0.5*(1/sqrt(pow(dx(x,y+1),2)+pow(dy(x,y+1),2)+pow(e,2))+ 1/sqrt(pow(dx(x,y),2)+pow(dy(x,y),2)+pow(e,2)));

			double g_4=0.5*(1/sqrt(pow(dx(x,y-1),2)+pow(dy(x,y-1),2)+pow(e,2))+ 1/sqrt(pow(dx(x,y),2)+pow(dy(x,y),2)+pow(e,2)));

			result(x,y)=(1-t*((g_1+g_2+g_3+g_4)))*aImage(x,y)+t*(g_1*aImage(x+1,y)+g_2*aImage(x-1,y)+g_3*aImage(x,y+1)+g_4*aImage(x,y-1));

		}
		aImage=result;

		if(k>1)
		{aImage=TVnlindiffusion(aImage, t, e,border,--k);}

		return aImage;


}

//nonlinear diffusion Potts potential
CMatrix<float> POTTSnlindiffusion(CMatrix<float> aImage, float t, int border, int k)
{

	CDerivative<float> aDerivative(3);
	CMatrix<float> dx(aImage.xSize(),aImage.ySize());
	CMatrix<float> dy(aImage.xSize(),aImage.ySize());
	NFilter::filter(aImage,dx,aDerivative,1);
	NFilter::filter(aImage,dy,1,aDerivative);
	CMatrix<float> result(aImage.xSize(),aImage.ySize(),255);


	for(int x=border;x<aImage.xSize()-border;x++)
		for(int y=border;y<aImage.ySize()-border;y++)
		{

			double g_1=0.5*( pow(dx(x+1,y),2)+pow(dy(x+1,y),2)+ pow(dx(x,y),2)+pow(dy(x,y),2));
			if(g_1==0)
			{
				g_1=0;
			}
			else
			{g_1=1;}

			double g_2=0.5*( pow(dx(x-1,y),2)+pow(dy(x-1,y),2)+ pow(dx(x,y),2)+pow(dy(x,y),2));
			if(g_2==0)
			{
				g_2=0;
			}
			else
			{g_2=1;}

			double g_3=0.5*( pow(dx(x,y+1),2)+pow(dy(x,y+1),2)+ pow(dx(x,y),2)+pow(dy(x,y),2));
			if(g_3==0)
			{
				g_3=0;
			}
			else
			{g_3=1;}

			double g_4=0.5*( pow(dx(x,y-1),2)+pow(dy(x,y-1),2)+ pow(dx(x,y),2)+pow(dy(x,y),2));
			if(g_4==0)
			{
				g_4=0;
			}
			else
			{g_4=1;}

			result(x,y)=(1-t*((g_1+g_2+g_3+g_4)))*aImage(x,y)+t*(g_1*aImage(x+1,y)+g_2*aImage(x-1,y)+g_3*aImage(x,y+1)+g_4*aImage(x,y-1));

		}
		aImage=result;

		if(k>1)
		{aImage=POTTSnlindiffusion(aImage, t,border,--k);}

		return aImage;


}

// nonlinear diffusion Perona-Malik

CMatrix<float> PMnlindiffusion(CMatrix<float> aImage, float t,float lambda, int border, int k){
	if(lambda==0)
	{
		lambda=0.01;
	}
	else
	{lambda=lambda;}
	CDerivative<float> aDerivative(3);
	CMatrix<float> dx(aImage.xSize(),aImage.ySize());
	CMatrix<float> dy(aImage.xSize(),aImage.ySize());
	NFilter::filter(aImage,dx,aDerivative,1);
	NFilter::filter(aImage,dy,1,aDerivative);
	CMatrix<float> result(aImage.xSize(),aImage.ySize(),255);


	for(int x=border;x<aImage.xSize()-border;x++)
		for(int y=border;y<aImage.ySize()-border;y++)
		{

			double g_1=0.5*(exp(-(pow(dx(x+1,y),2)+pow(dy(x+1,y),2))/pow(lambda,2))+ exp(-(pow(dx(x,y),2)+pow(dy(x,y),2))/pow(lambda,2)));

			double g_2=0.5*(exp(-(pow(dx(x-1,y),2)+pow(dy(x-1,y),2))/pow(lambda,2))+ exp(-(pow(dx(x,y),2)+pow(dy(x,y),2))/pow(lambda,2)));

			double g_3=0.5*(exp(-(pow(dx(x,y+1),2)+pow(dy(x,y+1),2))/pow(lambda,2))+ exp(-(pow(dx(x,y),2)+pow(dy(x,y),2))/pow(lambda,2)));

			double g_4=0.5*(exp(-(pow(dx(x,y-1),2)+pow(dy(x,y-1),2))/pow(lambda,2))+ exp(-(pow(dx(x,y),2)+pow(dy(x,y),2))/pow(lambda,2)));

			result(x,y)=(1-t*((g_1+g_2+g_3+g_4)))*aImage(x,y)+t*(g_1*aImage(x+1,y)+g_2*aImage(x-1,y)+g_3*aImage(x,y+1)+g_4*aImage(x,y-1));
		}
		aImage=result;

		if(k>1)
		{aImage=PMnlindiffusion(aImage, t, lambda,border,--k);}

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


void applyTVflow( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, float t, float e, int border, int k){

	//Neumann condition
	red= Neumann_bound(red,border);
	green= Neumann_bound(green,border);
	blue= Neumann_bound(blue,border);
	//apply diffusion
	red=TVnlindiffusion(red,t,e,border,k);
	green=TVnlindiffusion(green,t,e,border,k);
	blue=TVnlindiffusion(blue,t,e,border,k);

	cut(red,border);
	cut(green,border);
	cut(blue,border);
}



void applyPotts( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue,float t, int border, int k){



	//Neumann condition
	red= Neumann_bound(red,border);
	green= Neumann_bound(green,border);
	blue= Neumann_bound(blue,border);

	//apply PMnlindiffusion
	red=POTTSnlindiffusion(red,t, border, k);
	green=POTTSnlindiffusion(green,t, border, k);
	blue=POTTSnlindiffusion(blue, t, border, k);

	cut(red,border);
	cut(green,border);
	cut(blue,border);

}

void applyPMnlindiffusion( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, float t,float lambda, int border, int k){

	//Neumann condition
	red= Neumann_bound(red,border);
	green= Neumann_bound(green,border);
	blue= Neumann_bound(blue,border);
	//apply PMnlindiffusion
	red=PMnlindiffusion(red,t,lambda,border,k);
	green=PMnlindiffusion(green,t,lambda,border,k);
	blue=PMnlindiffusion(blue,t,lambda,border,k);

	cut(red,border);
	cut(green,border);
	cut(blue,border);

}                                            

void applyDual( CMatrix<float>& red, CMatrix<float>& green, CMatrix<float>& blue, float t,float alpha, int border, int k){

	//Neumann condition
	red= Neumann_bound(red,border);
	green= Neumann_bound(green,border);
	blue= Neumann_bound(blue,border);
	//apply PMnlindiffusion
	red=DUALnlindiffusion(red,t,alpha,border,k);
	green=DUALnlindiffusion(green,t,alpha,border,k);
	blue=DUALnlindiffusion(blue,t,alpha,border,k);

	cut(red,border);
	cut(green,border);
	cut(blue,border);
}

int main(int argc, char** argv){


	std::string fileNameInput;
	//parameters
	float t=0.25;
	float lambda=1;
	float e=1;
	float alpha=20;   
	int border=2;
	int k=20;
	int n;


	if (argc==2){
		fileNameInput=argv[1];

	} else if (argc==3){
		fileNameInput=argv[1];
		k=atoi(argv[2]);

	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		std::cout<<"Usage: gaussf inputfile <number of iteration k>"<<"\n";
		std::cout<<"The command should contain at least input file name. The default k=20."<<"\n";
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

	std::cout<<"Choose diffusion algorithm:"<<"\n";
	std::cout<<"[1]- TVflow; [2]- Dual; [3]- Perona-Malik; [4]-Potts potential"<<"\n";
	std::cin>>n;    
	if(n==1){
		applyTVflow(red, green, blue, t, e, border, k);
	}else if(n==2){
		applyDual(red,green, blue, t, alpha,  border, k);
	}else if(n==3){
		applyPMnlindiffusion( red, green,  blue,  t, lambda, border,  k);
	}else if(n==4){
		applyPotts( red,  green,  blue, t,  border, k);
	}else{
		std::cout<<"!!!WRONG INPUT!!!"<<"\n";
		return 0;
	}

	rgb2image(result_color,red, green, blue);

	result_color.writeToPPM((fileNameInput+"_diffusion.ppm").c_str());


	return 0;
}

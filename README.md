Computer vision: filters
====================================================
All filters can process only files in ppm format.

[![Build Status](https://travis-ci.org/Dtananaev/cv_filters.svg?branch=master)](https://travis-ci.org/Dtananaev/cv_filters)
[![BSD2 License](http://img.shields.io/badge/license-BSD2-brightgreen.svg)](https://github.com/Dtananaev/cv_filters/blob/master/LICENSE.md) 
     
To install type in the terminal in Linux:
* cd cv_filters
* mkdir build
* cd build
* cmake  ..
* make

It contains:

* gauss_noise - add gauss noise to the image 
<p align="center">
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/lamba.jpg" width="350"/>
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/lamborghini_noisy.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/gauss_noise
          * make
     * To run: ./gaussn name_of_file(without .ppm) \<mean value\> \<variance\>
     
* gauss_filter - add gauss blur to the image
<p align="center">
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla.jpg" width="350"/>
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla_gaussn_gaussf.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/gauss_filter
          * make
     * To run: ./gaussf inputfile(without .ppm) \<sigma of Gauss kernell\>
     
* box_filter - add blur with box filter
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/Chevrolet-Volt.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/volt_gaussn_boxf.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/box_filter
          * make
     * To run: ./boxf inputfile(without .ppm) \<size of box\> \<number of iterations of  box filter\>
     
* wavelet_decomposition - wavelet decomposition of the picture
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla_coeff.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/wavelet_decompose
          * make
     * To run: ./wavelet inputfile(without .ppm)  \<level of decomposition\> 
     
* wavelet_filter - wavelet schrinkage of the picture
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla_gaussn_coeff.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla_gaussn_recovered.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/wavelet_filter
          * make
     * To run: ./waveletf  inputfile(without .ppm)  \<level of decomposition\>  \<treshold for schrinkage\>
     * Additional information: it is possible to use 1 -hard schrinkage; 2 - soft schrinkage; 3 - Garrote schrinkage
     
* diffusion -  diffusion of the picture
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/Eiffel_tover.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tower_diffusion.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/diffusion
          * make
     * To run: ./diffusion  inputfile(without .ppm)  \<number of iteration k\>
     * Additional information: it is possible to use 1 -TVflow; 2 - Dual implementation; 3 -Perona-Malik diffusion; 4 - Potts potential

* nlm -  non-local mean filter
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/fallingMangoes_gaussn.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/fallingMangoes_gaussn_nlm.jpg" width="350"/>
</p>

     * To install use in terminal: 
          * cd ../cv_filter/nlm
          * make
     * To run: ./nlm  inputfile(without .ppm)  \<path radius\> \<window radius\> \<sigma\>

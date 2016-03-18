Computer vision: filters
====================================================

[![BSD2 License](http://img.shields.io/badge/license-BSD2-brightgreen.svg)](https://github.com/Dtananaev/cv_filters/blob/master/LICENSE.md) 

It contains:

* gauss_noise - add gauss noise to the image 
<p align="center">
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/lamba.jpg" width="350"/>
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/lamborghini_noisy.jpg" width="350"/>
</p>
     * To install use in terminal: 
          * cd ../cv_filter/gauss_noise
          * make
     * To run: ./gaussn name_of_file \<mean value\> \<variance\>
* gauss_filter - add gauss blur to the image
<p align="center">
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla.jpg" width="350"/>
  <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/tesla_gaussn_gaussf.jpg" width="350"/>
</p>
     * To install use in terminal: 
          * cd ../cv_filter/gauss_filter
          * make
     * To run: ./gaussf inputfile \<sigma of Gauss kernell\>
* box_filter - add blur with box filter
<p align="center">
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/Chevrolet-Volt.jpg" width="350"/>
   <img src="https://github.com/Dtananaev/cv_filters/blob/master/pictures/volt_gaussn_boxf.jpg" width="350"/>
</p>
     * To install use in terminal: 
          * cd ../cv_filter/box_filter
          * make
     * To run: ./boxf inputfile <size of box> <number of iterations of  box filter>

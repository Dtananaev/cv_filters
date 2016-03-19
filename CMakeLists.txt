cmake_minimum_required(VERSION 2.8)

project(cv_filters)

#Box filter
add_executable (boxf box_filter/boxf.cpp )
target_link_libraries (boxf)

#Gauss noise
add_executable (gaussn gauss_noise/gaussn.cpp )
target_link_libraries (gaussn)

#Gauss filter
add_executable (gaussf gauss_filter/gaussf.cpp )
target_link_libraries (gaussf)

#Diffusion filter
add_executable (diffusion diffusion/diffusion.cpp )
target_link_libraries (diffusion)

#Wavelet decompose
add_executable (wavelet wavelet_decompose/wavelet.cpp )
target_link_libraries (wavelet)

#Wavelet filter
add_executable (waveletf wavelet_filter/waveletf.cpp )
target_link_libraries (waveletf)

#Nlm filter
add_executable (nlm nlm_filter/nlm.cpp )
target_link_libraries (nlm)
 

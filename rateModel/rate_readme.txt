Contributed by Daniel Jercog.

This Matlab code generates r_E(t) and r_I(t) population rates and a adaptation a(t). Uses a mex-function eiModel (C++ source included, linux built mexa64).

___________________________

To compile the mex-function (tested on Matlab 2009 - Linux):
 >> mex -lgsl -lgslcblas eiModel.cpp

Example call to generate traces like those from Fig 5C:
 >> singleSimEIModel(0.7,3.5,4.8,20)


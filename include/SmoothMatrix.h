//================================================================
//
// Created: May 2001
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description: Thin plate spline energy functional
//                                                                           
//================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//================================================================
#ifndef _SMOOTHMATRIX_H
#define _SMOOTHMATRIX_H

/**
Note that the smooth matrix does not occupy much memory. Thus
it can be reinitialized many times during relaxion cycles.
*/
class SmoothMatrix {
  int K_; // order
  int n1_, n2_; // ??? NOTE: This is not the matrix dimension, but
  double alphaU_; // "alpha" , "dilation" factor (1/scale)
  double alphaU3_; // = alphaU_^3
  double alphaV_;
  double alphaV3_; // = alphaV_^3
  static double innerProduct_[3][4]; // [der][trans]

  static double innerProductQd0_[3][3]; // [offset-1][trans]
  static double innerProductQd1_[3][3]; // [offset-1][trans]
  static double innerProductQd2_[3][3]; // [offset-1][trans]
  int cutOff(int ij, int rs, int mn) const;

public:
  SmoothMatrix(int m, int n, double dU, double dV);
  ~SmoothMatrix(){}

  // indices runs from 0
  double operator()(int i, int j, int r, int s) const;
  double norm_l2() const;
};

#endif

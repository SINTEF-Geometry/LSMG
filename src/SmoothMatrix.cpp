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

#include <SmoothMatrix.h>

#ifdef WIN32
#include <minmax.h>
#endif

//#include <cmath>
#include <math.h>
#include <algorithm>
//#include <cstdlib>
#include <stdlib.h>

#include <iostream>
using namespace std;

// See also external documentation
double SmoothMatrix::innerProduct_[3][4] = {
  {151./315., 397./1680., 1./42., 1./5040.},
  {2./3., -1./8., -1./5., -1./120.},
  {8./3., -3./2., 0., 1./6.}
};

// Inner product betwen bases when they are cut off
// at the boundary; see external documentation,
// report by Øyvind Hjelle.
// The -99999999999999. are not used
double SmoothMatrix::innerProductQd0_[3][3] = {
  {599./1260.,-99999999999999.,-99999999999999999.},
  {151./630.,59./280.,-999999999999999.},
  {1./252.,43./1680.,1./84.}
};
double SmoothMatrix::innerProductQd1_[3][3] = {
  {37./60.,99999999999999.,99999999999999999.},
  {1./3.,-11./60.,999999999999999.},
  {1./20.,7./120.,-1./10.}
};
double SmoothMatrix::innerProductQd2_[3][3] = {
  {7./3.,99999999999999.,99999999999999999.},
  {4./3.,-1.,999999999999999.},
  {1./3.,-1./2.,0.}
};

//-----------------------------------------------------------------------------
SmoothMatrix::SmoothMatrix(int n1, int n2, double dU, double dV) {
  K_ = 4;  // order
  n1_ = n1;
  n2_ = n2;
  alphaU_ = 1./dU; // scaling/dilation in u-direction
  alphaU3_ = alphaU_*alphaU_*alphaU_;

  alphaV_ = 1./dV; // scaling/dilation in v-direction
  alphaV3_ = alphaV_*alphaV_*alphaV_;
}

//-----------------------------------------------------------------------------
int SmoothMatrix::cutOff(int ij, int rs, int n1n2) const {
  
  // ij can be i or j, ...., n1n2 can be n1 or n2

  // ????? Note also, if the grid is too small we have cut-off in both ends !!!!! 
  // ????? But this is not handled yet.
  // ????? So, assumes that m,n >= 3 such that cut off is only active in one end
  // See extended documentation for details

#ifdef DEBUG
  if (n1_ < 6 || n2_ < 6) { // ?????
    cout << "FATAL ERROR, too small grid for cut-off version" << endl;
    exit(-1);
  }
#endif

  int q = 3 - min(ij,rs);
  if (q > 0)
    return q;

  q = max(ij,rs) - (n1n2-4);
  if (q > 0)
    return q;

  return 0;
}

//-----------------------------------------------------------------------------
double SmoothMatrix::operator()(int i, int j, int r, int s) const {

  // Indices from 0

  int transU = abs(i-r);
  int transV = abs(j-s);

#if DEBUG
  if (transU >= K_ || transV >= K_) {
    //????
    cout << "WARNING: ZERO ELEMENT OF SmoothMatrix (i,j,r,s)= " << i << "," << j << "," << r << "," << s << endl;
    exit(-1);
    return 0.0;
  }
#endif
    
  // Find if we have a special case at the boundary in which case the bases are "cut off"
  // at the boundary of the domain:
  // if (i < 3 || j < 3 || r < 3 || s < 3 || i > n1-4 || j > n2-4 || r > n1-4 || s > n2-4)

  // No. "of cut-off intervals" in each dimension is now 1,2 or 3 (0 for no special case)
  // Let k be the translate between bases and let q be the cut-off (no. of intervals).
  // We have a special case if and only if q > k
  
  // We only need to calculate the cut-off and then use the look-up tables.
  // See extended documentation.

  // No. of cut-off intervals, > 0 gives special case
  // Here we give indices and dimension as in MBA
  // (runs from -1 up to m+1 and n+1)
  int qU = cutOff(i, r, n1_); // In U-dimension
  int qV = cutOff(j, s, n2_); // In V-dimension
    
  double innerU_d0, innerV_d0; // zero derivative
  double innerU_d1, innerV_d1; // first derivative
  double innerU_d2, innerV_d2; // second derivative

  
  // ????? must check if scalings are also valid in cut-off cases

  if (qU <= transU) {
    innerU_d0 = innerProduct_[0][transU]/alphaU_;
    innerU_d1 = innerProduct_[1][transU]*alphaU_;
    innerU_d2 = innerProduct_[2][transU]*alphaU3_;
  }
  else {
    innerU_d0 = innerProductQd0_[qU-1][transU]/alphaU_;
    innerU_d1 = innerProductQd1_[qU-1][transU]*alphaU_;
    innerU_d2 = innerProductQd2_[qU-1][transU]*alphaU3_;
  }

  if (qV <= transV) {
    innerV_d0 = innerProduct_[0][transV]/alphaV_;  // zero derivative
    innerV_d1 = innerProduct_[1][transV]*alphaV_;  // first derivative
    innerV_d2 = innerProduct_[2][transV]*alphaV3_; // second derivative
  }
  else {
    innerV_d0 = innerProductQd0_[qV-1][transV]/alphaV_;
    innerV_d1 = innerProductQd1_[qV-1][transV]*alphaV_;
    innerV_d2 = innerProductQd2_[qV-1][transV]*alphaV3_;
  }

  // return A + 2B + C
  // ???? Optimize this, e.g., by including the scalings in the matrix
  // But we can probably not include both alphaU_ and alphaV_ except for innerUandV_d1.
  // Anyway, we can save many multiplications...
  // - 2.0* can be included

  return innerU_d2 * innerV_d0 + 2.0*innerU_d1 * innerV_d1 + innerU_d0 * innerV_d2;

}

double SmoothMatrix::norm_l2() const {
  double sum2 = 0.0;
  int i,j,r,s;
  int L_ = K_;

  for(j = 0; j < n2_; j++) {
    int jmin = max(0,j-L_+1);
    int jmax = min(n2_-1,j+L_-1);
    for(i = 0; i < n1_; i++) {
      //int ii = j * n1_ + i; // row number
      int imin = max(0,i-K_+1);
      int imax = min(n1_-1,i+K_-1);
      for(s = jmin; s <= jmax; s++) {
        for(r = imin; r <= imax; r++) {
          double t = operator()(i,j,r,s);
          sum2 += t*t;
        }
      }
    }
  }
  return sqrt(sum2);
}

/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of LSMG.
 *
 * LSMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * LSMG is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with LSMG. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using LSMG.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the LSMG library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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

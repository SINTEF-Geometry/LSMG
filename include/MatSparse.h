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

#ifndef _MATSPARSE_H
#define _MATSPARSE_H


#include <BitMatrix.h>

#include <vector>
#include <fstream>

// ???? #define LSMG_real float

class MatSparse {
private:
  int noR_, noC_;         // matrix dimensions
  int p_;                 // number of non-zeros
  std::vector<int> irow_; // the offset in a_ of the first non-zero for each row.
                          // Thus, irow_(0) == 0 always.
                          // (extended by one; see constructor)
  std::vector<int> jcol_; // The column no. of a given offset
  std::vector<double> a_; // the p_ non-zero elements
public:
  MatSparse() : noR_(0), noC_(0), p_(0) {}

  MatSparse(int n1, int n2, int num_nonzero)
    : noR_(n1), noC_(n2), p_(num_nonzero), irow_(n1+1), jcol_(num_nonzero), a_(num_nonzero) {
    irow_[noR_] = num_nonzero;
  }
  
  ~MatSparse(){}

  void init(const BitMatrix& pattern); // Init the sparse structure from a bitmap pattern

  // Init the sparse structure from an entry-array, possibly sorted lexicographically
  void init(std::vector<std::pair<int,int> >& entries, int m, int n, bool sorted = false);

  void init(const MatSparse& mat, bool copyElements=false); // Serves as a "copy constructor" if copyElements=true

  void operator += (const MatSparse& mat); // Add a matrix (assumes sparsity pattern of mat is exactly like this).
  void operator *= (double fac);

  int noRows() const {return noR_;}
  int noColumns() const {return noC_;};
  
        double& operator () (int i, int j);       // inefficient !!!
  const double& operator () (int i, int j) const; // inefficient !!!

  double getVal(int i, int j) const; // as aove but returns zero when not in sparsity pattern inefficient !!!

  int getOffset(int i, int j) const; // and thus for access operator below

  void resize(int n1, int n2, int num_nonzero);

  int& irow(int k) {return irow_[k];} // 0 <= k <= noR_ index of a_ (last element is past-the-end)
  const int& irow(int k) const {return irow_[k];}
        int& jcol(int k) {return jcol_[k];} // 0 <= k <= p_-1
  const int& jcol(int k) const {return jcol_[k];}
        double& operator () (int k) {return a_[k];} // 0 <= k <= p_-1
  const double& operator () (int k) const {return a_[k];}

  double norml2() const;
  bool isSymmetric() const; // assumes m=n
  bool isDiagonallyDominant() const; 
  void printPattern(std::ostream& os) const;
  int noNonZeros() const {return p_;}

  void print(std::ostream& os) const; // debug
  void printMatlab(std::ostream& os) const; // debug
  void printPatternGnuplot(std::ostream& os) const; // debug
};
#endif

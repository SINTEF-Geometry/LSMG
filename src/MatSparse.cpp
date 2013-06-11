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

#include <MatSparse.h>
//#include <cmath>
#include <math.h>
#include <algorithm>
#include <iostream>


//-----------------------------------------------------------------------------
void MatSparse::init(std::vector<std::pair<int,int> >& entries, int m, int n, bool sorted) {

  // Init the sparse structure from an entry-array, possibly sorted lexicographically

  if (!sorted)
    std::sort(entries.begin(), entries.end());


  resize(m, n, entries.size());
  
  // ??? if the whole row is zero?
  std::fill(irow_.begin(), irow_.end(), -1);
  
  // Assume now that the array is lex-sorted
  int offset=0;
  int prev_row = -999;
  std::vector<std::pair<int,int> >::const_iterator it;
  for (it = entries.begin(); it != entries.end(); ++it) {
    int row = it->first;
    int col = it->second;
    if (row != prev_row) {
      irow_[row] = offset;
      prev_row = row;
    }
    jcol_[offset] = col;
    offset++;

  }
  irow_[noR_] = offset;
  
}

//-----------------------------------------------------------------------------
void MatSparse::init(const BitMatrix& pattern) {
  

  std::cout << "MatSparse::init(const BitMatrix& pattern) ..., but obsolete ??? " << std::endl; 


  // Init the sparse structure
  
  int m,n;
  pattern.size(m,n);

  // find no. of non-zeros in the matrix
  int nonzeros=0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (pattern.get(i, j))
	nonzeros++;

  //cout << "Here... " << m << ' ' << n  << endl;
  resize(m, n, nonzeros);


  int offset = 0;
  for (int i = 0; i < m; i++) { // each row
    irow_[i] = offset; // but if the whole row is zero?
    
    for (int j = 0; j < n; j++) { // examine all columns in current row
      if (pattern.get(i,j)) {
	jcol_[offset] = j;
	a_[offset] = 0.0; // The value
	offset++;
      }
    }
  }

  irow_[m] = offset;  // the last which is "passed the end"
}


//-----------------------------------------------------------------------------
void MatSparse::init(const MatSparse& mat, bool copyElements) { // Serves as a "copy constructor"
  
  noR_ = mat.noR_;
  noC_ = mat.noC_;
  p_   = mat.p_;
  irow_.resize(mat.irow_.size());
  std::copy(mat.irow_.begin(), mat.irow_.end(), irow_.begin()); 

  jcol_.resize(mat.jcol_.size());
  std::copy(mat.jcol_.begin(), mat.jcol_.end(), jcol_.begin()); 

  a_.resize(mat.a_.size());

  if (copyElements)
    std::copy(mat.a_.begin(), mat.a_.end(), a_.begin());
  else
    std::fill(a_.begin(), a_.end(), 0.0);
}


//-----------------------------------------------------------------------------
void MatSparse::operator += (const MatSparse& mat) { // Add a matrix (assumes sparsity pattern of mat is exactly like this).
  
  std::vector<double>::iterator it1;
  std::vector<double>::const_iterator it2;
  for (it1 = a_.begin(), it2 = mat.a_.begin(); it1 != a_.end(); ++it1, ++it2)
    *it1 += *it2;
}

void MatSparse::operator *= (double fac) {
  std::vector<double>::iterator it;
  for (it = a_.begin(); it != a_.end(); ++it)
    *it *= fac;
}


//-----------------------------------------------------------------------------
void MatSparse::resize(int n1, int n2, int num_nonzero) {
  noR_ = n1;
  noC_ = n2;
  p_ = num_nonzero;
  irow_.resize(n1+1);
  jcol_.resize(num_nonzero);
  a_.resize(num_nonzero);
  irow_[noR_] = num_nonzero;
}
 

/*
//-----------------------------------------------------------------------------
void MatSparse::prod(const PrVec& x, PrVec& y) const // Find y = Ax {

  int i,k;
  for(i=0; i<noR_; i++) {
    y(i) = 0.0;
    for(k = irow(i); k < irow(i+1); k++) {
      y(i) += (*this)(k) * x(jcol(k));
    }
  }
}
*/


//-----------------------------------------------------------------------------
double& MatSparse::operator () (int i, int j) {

  for(int k=irow(i); k<irow(i+1); k++) {
    if(jcol(k) == j)
      return a_[k];
  }

  std::cerr << "ERROR op1: (" << i << ',' << j << ") not in sparsity pattern" << std::endl;
  exit(-1);
}

//-----------------------------------------------------------------------------
const double& MatSparse::operator () (int i, int j) const {

  for(int k=irow(i); k<irow(i+1); k++) {
    if(jcol(k) == j)
      return a_[k];
  }
  std::cerr << "ERROR op2: (" << i << ',' << j << ") not in sparsity pattern" << std::endl;
  exit(-1);
}

//-----------------------------------------------------------------------------
double MatSparse::getVal(int i, int j) const { // inefficient
  // Same as above but returns zero if not in sparsity pattern
  for(int k=irow(i); k<irow(i+1); k++) {
    if(jcol(k) == j)
      return a_[k];
  }
  return 0.0;
}

//-----------------------------------------------------------------------------
int MatSparse::getOffset(int i, int j) const {
  for(int k=irow(i); k<irow(i+1); k++) {
    if(jcol(k) == j)
      return k;
  }
  std::cerr << "ERROR op2: (" << i << ',' << j << ") not in sparsity pattern" << std::endl;
  exit(-1);
}

//-----------------------------------------------------------------------------
double MatSparse::norml2() const {
  
  double sum2 = 0.0;
  for (int i = 0; i < noR_; i++) {
    // row by row
    for (int offset = irow(i); offset < irow(i+1); offset++)
      sum2 += operator()(offset)*operator()(offset);    
  }
  return sqrt(sum2);
}

// very innefficient, also tests all off-diagonal elements twice
//????
//-----------------------------------------------------------------------------
bool MatSparse::isSymmetric() const {

  if (noR_ != noC_) // if not quadratic
    return false;

  for (int ii = 0; ii < noR_; ii++) {
    for (int offset = irow(ii); offset < irow(ii+1); offset++) {
      int jj = jcol(offset);
      if ( operator()(offset) != operator()(jj,ii))
        return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
bool MatSparse::isDiagonallyDominant() const {

  for (int ii = 0; ii < noR_; ii++) {

    double sum = 0.0;
    double a_ii = 0.0; // since the diagonal element may be one of the zeros
    
    for (int offset = irow(ii); offset < irow(ii+1); offset++) {
      // Are we on the diagonal?
      int jj = jcol(offset); // = column no. of this offset
      if (jj == ii)
        a_ii = fabs(operator()(offset));
      else
        sum += fabs(operator()(offset));
    }
    if (sum > a_ii)
      return false;
  }
  return true;
}

//-----------------------------------------------------------------------------
void MatSparse::printPattern(std::ostream& os) const {
  for (int ii = 0; ii < noR_; ii++) {
    int offset = irow(ii);
    int jj = jcol(offset);
    for (int j = 0; j < noC_; j++) {
      if (j == jj) {
        os << '1';
        offset++;
        jj = jcol(offset);
      }
      else
        os << '0';
    }
    os << '\n';
  }
}

//-----------------------------------------------------------------------------
void MatSparse::print(std::ostream& os) const {
  for (int ii = 0; ii < noR_; ii++) {
    int offset = irow(ii);
    int jj = jcol(offset);
    for (int j = 0; j < noC_; j++) {
      if (j == jj) {
        os << a_[offset] << ' ';
        //os << '1';
        offset++;
        jj = jcol(offset);
      }
      else
        os << "XXX ";
    }
    os << '\n';
  }
}

//-----------------------------------------------------------------------------
void MatSparse::printMatlab(std::ostream& os) const {
  std::cout << "MatSparse: Printing file to Matlab format..." << std::endl;
  os << "A = [\n";
  for (int ii = 0; ii < noR_; ii++) {
    int offset = irow(ii);
    int jj = jcol(offset);
    for (int j = 0; j < noC_; j++) {
      if (j == jj) {
        os << a_[offset] << ' ';
        offset++;
        jj = jcol(offset);
      }
      else
        os << "0 ";
    }
    os << ";\n";
  }
  os << ']';
}

//-----------------------------------------------------------------------------
void MatSparse::printPatternGnuplot(std::ostream& os) const {
  std::cout << "MatSparse: Printing file to Gnuplot format..." << std::endl;
  for (int ii = 0; ii < noR_; ii++) {
    int offset = irow(ii);
    int jj = jcol(offset);
    for (int j = 0; j < noC_; j++) {
      if (j == jj) {
        os << ii << ' ' << jj << '\n';
        offset++;
        jj = jcol(offset);
      }
      //else
      //os << "0 ";
    }
    //os << ";\n";
  }
}



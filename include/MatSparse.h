//================================================================
//
// Created: May 2001
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//         Originates from Michael Floater.
//                                                                           
// Revised:
//                                                                           
// Description: Representation of sparse matrix.
//              This version is not symmetric. 
//                                                                           
//================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//================================================================
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

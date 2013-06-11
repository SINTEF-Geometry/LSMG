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

#include <BitMatrix.h>

#include <iostream>
using namespace std;


//#define ARRAY_RANGE_CHECK 1


BitMatrix::BitMatrix(unsigned int sizx, unsigned int sizy) {
  matrix_=new unsigned char[(sizx/8+1)*sizy];
  sizex_=sizx;
  sizey_=sizy;
  zero();
}

BitMatrix::BitMatrix() {
  matrix_=0;
  sizex_=0;
  sizey_=0;
}

/*
BitMatrix::BitMatrix(BitMatrix &bm) {
	matrix=bm.matrix;
	sizex=bm.sizex;
	sizey=bm.sizey;
}
*/

BitMatrix::~BitMatrix() {
  if (matrix_!=0)
    delete matrix_;
}

void BitMatrix::set(const unsigned int& x, const unsigned int& y, int val) {

#ifdef ARRAY_RANGE_CHECK
  if (x<0 || x >=sizex_ || y<=0 || y>=sizey_) {
    cout << "BitMatrix::set(...), (x,y) outside range" << endl;
    exit(-1);
  }
#endif
  
  if (val==0)
    matrix_[y*(sizex_/8+1)+(x/8)]=matrix_[y*(sizex_/8+1)+(x/8)]&(~(1<<(x%8)));
  else
    matrix_[y*(sizex_/8+1)+(x/8)]=matrix_[y*(sizex_/8+1)+(x/8)]|(1<<(x%8));
}

int BitMatrix::get(const unsigned int& x, const unsigned int& y) const {
#ifdef ARRAY_RANGE_CHECK
  if (x<0 || x >=sizex_ || y<=0 || y>=sizey_) {
    cout << "BitMatrix::set(...), (x,y) outside range" << endl;
    exit(-1);
  }
#endif

  return ((matrix_[y*(sizex_/8+1)+(x/8)]&(1<<(x%8)))>>(x%8));
}

/*
void BitMatrix::raw() const {
  int c;
  for (c=0;c<(sizex_/8+1)*sizey_;c++)
    //printf("%X ",matrix_[c]);
  cout << '\n';
  return;
}
*/

void BitMatrix::zero() {
  int c;
  for (c=0;c<(sizex_/8+1)*sizey_;c++)
    matrix_[c]=0;
  return;
}

void BitMatrix::show() const {
  int x,y;
  for (y=0;y<sizey_;y++) {
    for (x=0;x<sizex_;x++) {
      //printf("%1i ",get(x,y));
      cout << get(x,y) << ' ';
    }
    cout << '\n';
  }
  return;
}

unsigned long BitMatrix::count(int val) const {
  int x,y;
  unsigned long int c=0;
  for (y=0;y<sizey_;y++) {
    for (x=0;x<sizex_;x++) {
      if (get(x,y)==val)
	c++;
    }
  }
  return c;
}

// -----------------------------------------------------------------------

//#define TEST_BITMATRIX
#ifdef TEST_BITMATRIX
void main() {
  int nx=10; int ny=10;
  BitMatrix mat(nx,ny);

  cout << "raw:" << endl;
  //mat.raw();

  mat.set(0,0,1);
  mat.set(2,0,1);
  mat.set(nx-1,ny-1,1);
  mat.set(4,7,1);
  cout << "show1:" << endl;
  mat.show();

  int val = 0; cout << "No. of vals: " << val << '=' << mat.count(val) << endl;
}
#endif

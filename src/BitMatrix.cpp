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

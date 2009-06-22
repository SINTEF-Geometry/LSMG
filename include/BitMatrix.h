#ifndef _BITMATRIX_H_
#define _BITMATRIX_H_


// Found the basics at Björn Arne Herwig's homepage:
// "Mein wahrscheinlich erstes C++-Programm"
class BitMatrix {
  
  unsigned int sizex_, sizey_;
  unsigned char* matrix_;
 public:
  BitMatrix(unsigned int sizx, unsigned int sizy);
  BitMatrix();
  ~BitMatrix();
 
  void set(const unsigned int& x, const unsigned int& y, int val);
  int  get(const unsigned int& x, const unsigned int& y) const;
  
  void size(int& nox, int& noy) const {nox=sizex_; noy=sizey_;} 

  //void raw() const;
  void zero();
  void show() const;
  unsigned long count(int val) const;
};

#endif

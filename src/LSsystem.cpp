//================================================================
//
// Created: September 2001
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description: Sparse solver
//                                                                           
//================================================================
// Copyright (c) 2002 SINTEF Applied Mathematics
//================================================================
#include <LSsystem.h>
#include <stdexcept>

#ifdef WIN32
#include <minmax.h>
#endif

#include <SmoothMatrix.h>

#include <UCBsplines.h>
#include <MBAclock.h> // ???

#define DEBUG_LSMG_0  1    // output such as number of iterations
//#define DEBUG_LSMG_1  1 // misc. extra checks and information

#include <algorithm>
using namespace std;

#include <MGsystem.h> // This is only needed for relaxCGPrecond where we use MG for preconditioning

static double average(std::vector<double> v) {
  
    std::vector<double>::const_iterator it;
    double sum = 0.0;
    for (it = v.begin(); it != v.end(); ++it) {
	sum += *it;
    }
    return sum/v.size();
}

//-------------------------------------------------------------------------
LSsystem::LSsystem(boost::shared_ptr<std::vector<double> > U,
                   boost::shared_ptr<std::vector<double> > V,
                   boost::shared_ptr<std::vector<double> > Zvals,
                   int noX, int noY, // Thus, number of unknowns is noX*noY
                   boost::shared_ptr<GenMatrix<UCBspl_real> > x,
                   double lambdaFac)
    : U_(U), V_(V), Zvals_(Zvals), x_(x) {

    n1_ = noX ; n2_ = noY; // This is also done when building equation system
    K_=L_=4;
    lambdaFac_ = lambdaFac;

    if (!x_.get())
	x_.reset(new GenMatrix<UCBspl_real> );

    if (x_->noX() != noX || x_->noY() != noY) {
	x_->resize(noX,noY);
	x_->fill(average(*Zvals));
    }
}


//-------------------------------------------------------------------------
void LSsystem::calculateDomain() {
  
    if (domain_.size() != 4)
	domain_.resize(4);

    domain_[0] = *std::min_element(U_->begin(), U_->end());
    domain_[1] = *std::min_element(V_->begin(), V_->end());
    domain_[2] = *std::max_element(U_->begin(), U_->end());
    domain_[3] = *std::max_element(V_->begin(), V_->end());
}

//-------------------------------------------------------------------------
void LSsystem::setDomain(double umin, double vmin, double umax, double vmax) {
    if (domain_.size() != 4)
	domain_.resize(4);

    domain_[0] = umin;
    domain_[1] = vmin;
    domain_[2] = umax;
    domain_[3] = vmax;
}

//-------------------------------------------------------------------------
void LSsystem::getDomain(double& umin, double& vmin, double& umax, double& vmax) {

    umin = domain_[0];
    vmin = domain_[1];
    umax = domain_[2];
    vmax = domain_[3];
}


//-------------------------------------------------------------------------
static int noNonZeros(int n1_, int n2_, int K_, int L_) {
    // Count non-zeros
    int numNonZeros = 0;
    int imin,imax,jmin,jmax;
    for(int j = 0; j < n2_; j++) {
	jmin = max(0,j-L_+1);
	jmax = min(n2_-1,j+L_-1);
    
	for(int i = 0; i < n1_; i++) {
	    imin = max(0,i-K_+1);
	    imax = min(n1_-1,i+K_-1);
	    numNonZeros += (imax - imin + 1) * (jmax - jmin +1);
	}
    }
    return numNonZeros;
}

//-------------------------------------------------------------------------
void LSsystem::buildSparseStructure() {

    // Must be run after n1_ and n2_ have been set
  
    if (n1_ < 1 || n2_ < 1) {
	throw runtime_error("ERROR in buildSparseStructure, inproper initialization");
    }

    // Build the data structure for A_ and init with zeros
  
    // Count non-zeros
    int numNonZeros = noNonZeros(n1_, n2_, K_, L_);
  
    A_.resize(n1_*n2_, n1_*n2_, numNonZeros);
  
    // Fill in data structure for sparse matrix
    int offset = 0; // Index from zero
    for(int j = 0; j < n2_; j++) {
	int jmin = max(0,j-L_+1);
	int jmax = min(n2_-1,j+L_-1);
	for(int i = 0; i < n1_; i++) {
	    int ii = j * n1_ + i; // row number
	    A_.irow(ii) = offset; // Note, now it goes from 0
	    int imin = max(0,i-K_+1);
	    int imax = min(n1_-1,i+K_-1);
	    for(int s = jmin; s <= jmax; s++) {
		for(int r = imin; r <= imax; r++) {
		    A_.jcol(offset) = s * n1_ + r;
		    A_(offset) = 0.0;
		    offset++;
		}
	    }
	}
    }
    A_.irow(n1_*n2_) = offset;
}

//------------------------------------------------------------------------
void LSsystem::calcLambda() {
  
    // Calculate lambda from the norm of the
    // matrices: lambda = Norm_A/norm_E where E is the smoothing matrix.
    // This value will be multiplied with lambdaFac_ (if lamdaFac_ is different
    // from 1.0, which is the default).
    // Must  be done after A_ is calculated (from the scattered data)
    // (Note that after addPoint, the smoothing term is not correct.)
  
    // Calculate lambda for smoothing term
    double sa_norm = A_.norml2();
  
    int m_ = n1_-3;
    int n_ = n2_-3;
  
    double du = (umax() - umin()) / (double)m_;
    double dv = (vmax() - vmin()) / (double)n_;
  
    SmoothMatrix E(m_+3,n_+3,du,dv);
  
    double e_norm = E.norm_l2();
    lambda_ = 1.0;
    if (e_norm != 0.0) {
	lambda_ = sa_norm/e_norm;
    } else {
	throw runtime_error("LSsystem::calcLambda(): e_norm is zero when "
			    "calculating lambda, exit now in Debug mode");
    }
}

//------------------------------------------------------------------------
void LSsystem::addSmoothingTerm(double lambdaVal) {
  
    // where lambdaVal is either the original lambda_*lambdaFac_,
    // which is given when the system is built, or it is adjustment
    // if we want to smooth/recover the surface later.
  
    // Note that the smoothing term does not depend on the scattered data,
    // Only the system matrix is affected.
  
    int m_ = n1_-3;
    int n_ = n2_-3;  
    double du = (umax() - umin()) / (double)m_;
    double dv = (vmax() - vmin()) / (double)n_;
    SmoothMatrix E(m_+3,n_+3,du,dv);
  
    // ADD lambda * E to G.
    int noRows = n1_*n2_;
    for(int ii = 0; ii < noRows ; ii++) {
	int i = ii % n1_;
	int j = ii/n1_;
	for(int offset = A_.irow(ii); offset < A_.irow(ii+1); offset++) {
	    int jj = A_.jcol(offset); // column no. in sA_
	    int r = jj % n1_;
	    int s = jj/n1_;
	    A_(offset) += lambdaVal * E(i,j,r,s);
	}
    }
}

//------------------------------------------------------------------------
void LSsystem::setSmoothingFactor(double newLambdaFac) {
    
    // adjustment such that we can use the same matrix
    double lambdaAdjust = lambda_*(newLambdaFac - lambdaFac_);
  
    lambdaFac_ = newLambdaFac;

    addSmoothingTerm(lambdaAdjust);
}



//------------------------------------------------------------------------
void LSsystem::buildRHS() {
  
    // Build the right hand side (RHS) of the equation system for the input data,
    // including the z-values which are only used here.
    // The code is extracted from buildEqSystem below which also optionally
    // builds the RHS.  

    if (domain_.size() != 4)
	calculateDomain();

    // NOTE: it is assumed that the dimension is present in n1_ and n2_
    b_.resize(n1_, n2_);
    b_.fill(0.0);
  
    int noPoints = U_->size();  
    int m_ = n1_-3;
    int n_ = n2_-3;
    double u_scaling = double(m_) / (umax() - umin());
    double v_scaling = double(n_) / (vmax() - vmin());

    int i, j, k, l;
    double s, t;
    double w_kl[4][4];

    if (weights_.get() && weights_->size() != 0) {
	// weighed case
	for (int ip = 0; ip < noPoints; ip++) {
	    // Map to the half open domain Omega = [0,m) x [0,n)
	    // The mapped uc and vc must be (strictly) less than m and n respectively
	    double uc = ((*U_)[ip] - umin()) * u_scaling;
	    double vc = ((*V_)[ip] - vmin()) * v_scaling;
	    
	    UCBspl::ijst(m_, n_, uc, vc, i, j, s, t); // i and j from -1
	    
	    // All the 16 tensors in the 4x4 neighbourhood of a point
	    UCBspl::WKL(s, t, w_kl); // substituted by Odd Andersen, 16 dec. 2003
	    
	    double zw = (*Zvals_)[ip];
	    
	    map<int, double, less<int> >::iterator mit = weights_->find(ip);
	    if(mit != weights_->end()) 
		zw *= (*mit).second;
	    for (k = 0; k <= 3; k++) {
		for (l = 0; l <=3; l++) {
		    b_(i+k,j+l) += w_kl[k][l] * zw;
		}
	    }
	}
    } else {
	// unweighed case - remove weighed multiplication in loop
	// (this test is done outside the loop in order to speed up execution
	for (int ip = 0; ip < noPoints; ip++) {
	    // Map to the half open domain Omega = [0,m) x [0,n)
	    // The mapped uc and vc must be (strictly) less than m and n respectively
	    double uc = ((*U_)[ip] - umin()) * u_scaling;
	    double vc = ((*V_)[ip] - vmin()) * v_scaling;
	    
	    UCBspl::ijst(m_, n_, uc, vc, i, j, s, t); // i and j from -1
	    
	    // All the 16 tensors in the 4x4 neighbourhood of a point
	    UCBspl::WKL(s, t, w_kl); // substituted by Odd Andersen, 16 dec. 2003
	    
	    double zw = (*Zvals_)[ip];
	    
	    for (k = 0; k <= 3; k++) {
		for (l = 0; l <=3; l++) {
		    b_(i+k,j+l) += w_kl[k][l] * zw;
		}
	    }
	}
    }
}


//-----------------------------------------------------------------------------------
bool LSsystem::addPoint(double u, double v, double z, double weight, bool addToRHS) 
{  
    // (Observe: z is dummy if addToRHS is false)
  
    // Add point to system matrix and RHS
    int m_ = n1_-3;
    int n_ = n2_-3;
  
    // Map to the half open domain Omega = [0,m) x [0,n)
    // The mapped uc and vc must be (strictly) less than m and n respectively
    double uc = (u - umin())/(umax()-umin()) * (double)m_;
    double vc = (v - vmin())/(vmax()-vmin()) * (double)n_;
  
    if (uc < 0.0 || vc < 0.0 || uc > m_ || vc > n_) {
	throw runtime_error("ERROR in LSsystem::buildRHS: A point was mapped "
			    "to outside domain");
    }
  
    int i, j;
    double s, t;
    UCBspl::ijst(m_, n_, uc, vc, i, j, s, t); // i and j from -1
  
    double w_kl[4][4];
  
    // All the 16 tensors in the 4x4 neighbourhood of a point
    UCBspl::WKL(s, t, w_kl); // substituted by Odd Andersen, 16 dec. 2003

    int k, l;
    double weighted_z = z * weight;
    if (addToRHS) {
	// The right hand side (RHS) elements
	// ----------------------------------
	for (k = 0; k <= 3; k++) {
	    for (l = 0; l <=3; l++) {
		b_(i+k,j+l) += w_kl[k][l] * weighted_z;
	    }
	}
    }
    
    // The sparse system matrix
    for(k = 0; k < K_; k++) {
	for(l = 0; l < L_; l++) {
	    // This gives: w_kl[k][l]
      
	    // i and j are calculated above and goes from -1
	    // while the sparse matrix goes from 0:    
      
	    // row no. in sparse matrix (linearized 1D index from 2D index)
	    // ii runs from 0
	    int ii = (j + l + 1) * n1_ + (i + k + 1); // Note + L + 1
      
	    int imin = max(0,i+k - K_ + 2);
	    int jmin = max(0,j+l - L_ + 2);
	    int imax = min(n1_-1, i + k + K_ - 2 + 2);
      
	    // Then, the tensor to take the product with (also utilize symmetry?)
	    
	    for(int kk = 0; kk < K_; kk++) {
		for(int ll = 0; ll < L_; ll++) {
          
		    int r = i+kk - K_ + 5;
		    int s = j+ll - L_ + 5;
          
		    int offset = A_.irow(ii) + (s - jmin)*(imax - imin + 1) + (r - imin);
		    
		    double contrib = w_kl[k][l]*w_kl[kk][ll] * weight;
		    A_(offset) += contrib;
		}
	    }
	}
    } // end for k=0 ... 
  
    return true;
}


//------------------------------------------------------------------------
void LSsystem::buildEqSystem(bool buildRHS) {
  

    MBAclock rolex;

    if (domain_.size() != 4)
	calculateDomain();

    // In case the system has been updated by a friend (Multigrid)
    // so NOTE: it is assumed that the solution vector is present with correct dimensions
    n1_ = x_->noX() ; n2_ = x_->noY();

    // Build the data structure for the sparse matrix and fill it with zeros
    buildSparseStructure();
  
    // Fill in right hand side b_ and the sparse matrix A_
    // (The smoothing term will be added afterwards, see end of function)
    // Note that we also have a "buildRHS" above that can be used separately (in Multigrid)
    // ------------------------------------------------------------------------------------
  
    if (buildRHS) {
	b_.resize(n1_, n2_);
	b_.fill(0.0);
    }
  
    int noPoints = U_->size();  
    double weight = 1.0;
    bool wls = (weights_.get() && weights_->size() != 0);

    map<int,double>::iterator mit;
    for (int ip = 0; ip < noPoints; ip++) {
	if (wls) {

	    mit = weights_->find(ip);
	    if (mit != weights_->end())
		weight = (*mit).second;
	    else
		weight = 1.0;
	}      
	addPoint((*U_)[ip], (*V_)[ip], (*Zvals_)[ip], weight, buildRHS);
    }
  
    calcLambda(); // and set private member
    addSmoothingTerm(lambda_*lambdaFac_);
}

//------------------------------------------------------------------------
double LSsystem::rowMatVecMult(int row_no, const MatSparse& A, const GenMatrix<UCBspl_real>& x) {

    // STATIC class function (but can also be a free function)
  
    // Vector-vector product between the i'th row-vector in A with the vector x
    // (The vector x is represented as a matrix)
  
    // Matrix-vector product: Ax for row no. i
    double sum = 0.0;
  
    int n1 = x.noX();
    for (int offset = A.irow(row_no); offset < A.irow(row_no+1); offset++) { // for each row
    
	int j = A.jcol(offset); // the column number
	int k = j%n1 - 1; // runs fastest (-1,...,n1-2)
	int l = j/n1 - 1;
    
	sum += A(offset)*x(k,l);
    }
    return sum;
}

//------------------------------------------------------------------------
double LSsystem::currentSolutionNorm() const
{
    return sqrt((*x_).norm_2());
}

//------------------------------------------------------------------------
double LSsystem::residual_l2(bool scaled) const {
  
    // the l_2 norm of the residual; scaled if indicated by parameter scaled

    double h=1.0;
    if (scaled) {
	int m_ = n1_-3;
	int n_ = n2_-3;

	// scaling factor as the area of a grid cell
	double du = (umax() - umin()) / (double)m_;
	double dv = (vmax() - vmin()) / (double)n_;    
	h = du*dv;
    }  

    int n = n1_*n2_; // The "one-dimensional" size of RHS and the soluton vector
    double norm2 = 0.0;
  
    for (int i = 0; i < n; i++) { // for each row no. i in the system matrix
    
	// Matrix-vector product: Ax
	double sum = rowMatVecMult(i,A_,*x_);
    
	// Determine indices of b
	int k = i%n1_ - 1; // runs fastest (-1,...,n1-2)
	int l = i/n1_ - 1;
    
	double val = b_(k,l) - sum; // b - Ax  element

	norm2 += val * val;

    }
    
    if (scaled) {
	norm2 *= h;
    }
  
    return sqrt(norm2);
}

//------------------------------------------------------------------------
double LSsystem::residual_linf() const {
  
    // the max norm (Chebychev norm) of the residual, i.e.,
    //
    // ||b - Ax||_inf
  
    int n = n1_*n2_; // The "one-dimensional" size of RHS and the soluton vector
    double norm = 0.0;
  
    for (int i = 0; i < n; i++) { // for each row no. i in the system matrix
    
	// Matrix-vector product: Ax
	double sum = rowMatVecMult(i,A_,*x_);
    
	// Determine indices of b
	int k = i%n1_ - 1; // runs fastest (-1,...,n1-2)
	int l = i/n1_ - 1;
    
	double val = fabs(b_(k,l) - sum); // b - Ax  element
	if (val > norm)
	    norm = val;
    }
  
    return norm;
}


//------------------------------------------------------------------------
// It is implemented in the friend class MGsystem
void LSsystem::residual(GenMatrix<UCBspl_real>& r) const {
  
    // NOTE: This is almost the same function as above
    // The residual will be represented as a matrix.
  
    //
    // b - Ax
  
    if (r.noX() != n1_ || r.noY() != n2_)
	r.resize(n1_, n2_);
    r.fill(0.0); // ??? not necessary
  
    int n = n1_*n2_; // The "one-dimensional" size of RHS and soluton vector
  
    for (int i = 0; i < n; i++) { // for each row no. i in the system matrix
    
	// Matrix-vector product: Ax for row no. i
	double sum = rowMatVecMult(i,A_,*x_);
    
	// Determine indices of b
	int k = i%n1_ - 1; // runs fastest (-1,...,n1-2)
	int l = i/n1_ - 1;
    
	r(k,l) = b_(k,l) - sum; // b - Ax  element
    }
}


//------------------------------------------------------------------------
static double norm2(GenMatrix<UCBspl_real>& v) {
  
    // Norm squared
    int noX = v.noX();
    int noY = v.noY();

    double sp = 0.0;
    for (int j = -1; j <= noY-2; j++) {
	for (int i = -1; i <= noX-2; i++) {
	    double tmp = v(i,j);
	    sp += tmp*tmp;
	}
    }
    return sp;
}



//------------------------------------------------------------------------
static double scalarProduct(GenMatrix<UCBspl_real>& v1, GenMatrix<UCBspl_real>& v2) {

    if (&v1 == &v2)
	return norm2(v1);

    int noX = v1.noX();
    int noY = v1.noY();

    double sp = 0.0;
    for (int j = -1; j <= noY-2; j++) {
	for (int i = -1; i <= noX-2; i++) {
	    sp += v1(i,j)*v2(i,j);
	}
    }
    return sp;
}


//------------------------------------------------------------------------------------
static void matVecMult(const MatSparse& A, const GenMatrix<UCBspl_real>& P,
		       GenMatrix<UCBspl_real>& T) {

    // STATIC class function (but can also be a free function).
    // Uses the static class function LSsystem::rowMatVecMult
    
    // Matrix-vector product between A and P, where also P is represented as a matrix, and
    // so is the resulting vector T

    int n1 = P.noX();
    int n2 = P.noY();
    if (T.noX() != n1 || T.noY() != n2)
	T.resize(n1, n2);
 
    int j = 0;

    int dimA = A.noRows();
    for (int row_no = 0; row_no < dimA; row_no++) {
	int k = j%n1 - 1; // runs fastest (-1,...,n1-2); see also LSsystem::rowMatVec
	int l = j/n1 - 1;
	T(k,l) = LSsystem::rowMatVecMult(row_no, A, P);
	j++;
    }
}


//------------------------------------------------------------------------
static void vecPlusScalVec_1(GenMatrix<UCBspl_real>& G, GenMatrix<UCBspl_real>& H,
			     double scalar) {

    // Calculates G + aplha_k*H  with G OVERWRITTEN

    int noX = G.noX();
    int noY = G.noY();

    // ??? change later to run el by el if genmatrix is changed?
    for (int j = -1; j <= noY-2; j++) {
	for (int i = -1; i <= noX-2; i++) {
	    G(i,j) += scalar*H(i,j);
	}
    }
}
static void vecPlusScalVec_2(GenMatrix<UCBspl_real>& G, GenMatrix<UCBspl_real>& H,
			     double scalar) {

    // Calculates G + aplha_k*H with H OVERWRITTEN

    int noX = G.noX();
    int noY = G.noY();

    // ??? change later to run el by el if genmatrix is changed?
    for (int j = -1; j <= noY-2; j++) {
	for (int i = -1; i <= noX-2; i++) {
	    double tmp = G(i,j) + scalar*H(i,j);
	    H(i,j) = tmp;
	}
    }
}


//------------------------------------------------------------------------
void LSsystem::relaxGaussSeidel(int noIterations) 
{
    if (noIterations == 0)
	return;

    int j;
    int num_unknowns = numberOfUnknowns();

    GenMatrix<UCBspl_real>& x = *x_;

    for (int n_it = 1; n_it <= noIterations; n_it++) {
    
	for (int i = 0; i < num_unknowns; i++) { // for each row no. i
	    double sumL = 0.0;
      
	    int offset; 
	    for (offset = A_.irow(i); (j=A_.jcol(offset)) < i; offset++) {
		// Find row and column in x (PHI_) from row and column in A
		// (The size of x is  n1 x n2  while the size of A is (n1*n2) x (n1*n2)
        
		int l = j/n1_;// - 1;
		int k = j - (n1_ * l) - 1;  // k = j%n1_ - 1;
		l -= 1;
        
		sumL += A_(offset) * x(k,l);
	    }
      
	    double sumR = 0.0;
      
	    int offset_diag = offset;
	    offset++; // step over diagonal
      
	    for ( ;offset < A_.irow(i+1); offset++) {
		j = A_.jcol(offset);

		int l = j/n1_;// - 1;
		int k = j - (n1_ * l) - 1;  // k = j%n1_ - 1;
		l -= 1;

		sumR += A_(offset) * x(k,l);
	    }
      
	    // Diagonal element in x (and b) is now determined by row no. i in A
	    // Lexicographically or 

	    int l = i/n1_;
	    int k = i - (n1_ * l) - 1; //k = i%n1_ - 1; // runs fastest (-1,...,n1-2)
	    l -= 1;

	    x(k,l) = (b_(k,l) - sumL - sumR)/A_(offset_diag);
      
	    // SOR method
	    // BUT:
	    // Note that SOR should not be used as a smoothing operator in MG.
	    // The overrelaxation destroys the high-frequency smoothing that
	    // is so crucial for the multigrid method !!!
	    //double x1 = (b_(k,l) - sumL - sumR)/A_(offset_diag);
	    //double x2 = x_(k,l);
	    //const double omega = 2.0/3.0;
	    //x_(k,l) = omega*x1 + (1.0 - omega)*x2;
      
	} // for rows
    } // iteration
}

//------------------------------------------------------------------------
void LSsystem::relaxCG(int noIterations) {

    // Conjugate gradient method, relax on Ax = f
  
    // ??? is float sufficient now?

    // See lecture notes (extract), Chapter 3 by Tom Lyche

    // Step1: Starting vector (defines the first residual)
    //----------------------------------------------------

    // Step 2: p_0 = r_0 = f - Ax_0
    GenMatrix<UCBspl_real> r_k(n1_,n2_); // = p0
    residual(r_k);

    GenMatrix<UCBspl_real> p_k;
    p_k.init(r_k);

    // Step 3: rho_0 = r_0^T * r_0
    double rho_k = scalarProduct(r_k,r_k);

    int k=0;
  
    // Step 4 (loop, starting with k=0)
    // --------------------------------   
    GenMatrix<UCBspl_real> t_k(n1_,n2_);


    // To avoid division by zero, we break the loop if necessary by testing againts eps_num.
    // (e.g., when the starting vector corresponds to the exact solution and the initial residual is zero).

    double eps_num = 1.0e-15;
    while (k < noIterations && rho_k > eps_num) {

	// Step 4.1: t_k = A * p_k
	::matVecMult(A_, p_k, t_k);
  
	// Step 4.2: alpha_k = rho_k/(p_k^T*t_k)
	double alpha_k = rho_k/scalarProduct(p_k,t_k);

	// Step 4.3: x_(k+1) = x_k + alpha_k * p_k (x_k can be overwritten)
	// (x_k is the algorithm is the private solution vector x_)
	::vecPlusScalVec_1(*x_, p_k, alpha_k);

	// Step 4.4: r_(k+1) = r_k - alpha_k * t_k (r_k can be overwritten)
	::vecPlusScalVec_1(r_k, t_k, -1.0*alpha_k);

	// Step 4.5: rho_(k+1) = r_(k+1)^T * r_(k+1)
	double rho_kPlus1 = scalarProduct(r_k,r_k);


	// Step 4.6: p_(k+1) = r_(k+1) + rho_(k+1)/rho_k * p_k (p_k can be overwritten)
	double rho_fac = rho_kPlus1/rho_k;
	::vecPlusScalVec_2(r_k, p_k, rho_fac);

	// Step 4.7:
	k += 1;

	// Update to next iteration
	// x_k, r_k and p_k have been overwritten and are ok
	rho_k = rho_kPlus1;
    }
}


//------------------------------------------------------------------------
void LSsystem::relaxCGPrecond(int noIterations) {

    // Conjugate gradient method with preconditioning.
    // In general we relax on the system BAx = Bf.
    // Here B^{-1} \approx A. But B is not available here. There are
    // two statements in the algorithm involving B on the right hand side.
    // We move B to the left hand side and and set it equal to A and run a
    // multigrid cycle on the system. ?????????????????

    // Using Multigrid as a preconditioner.
  
    // ??? is float sufficient now?

    // ??? domain

    // See lecture notes (extract), Chapter 3 by Tom Lyche

    // Step1: Starting vector (defines the first residual)
    //----------------------------------------------------

    // Step 2: p_0 = r_0 = f - Ax_0
    GenMatrix<UCBspl_real> r_k(n1_,n2_); // = p0
    residual(r_k);
  
    // This is the residual in the preconditioned BAx = Bf
    boost::shared_ptr<GenMatrix<UCBspl_real> > s_k_temp(new GenMatrix<UCBspl_real> (n1_,n2_));

    // Run an MG cycle on the system As_k = r_k to oblain s_k

    //------------------------------
    s_k_temp->fill(0.0);
    MGsystem mgsys;
    boost::shared_ptr<std::vector<double> > dummy1;
    boost::shared_ptr<std::map<int,double> > dummy2;
    mgsys.initMultilevel(1,1,6,U_,V_,dummy1,s_k_temp,dummy2); // ??? hard coded h=? now 
    mgsys.qweSetRHS(r_k); // Set the initial right hand side at the finest level. ??? temporary function now
    //mgsys.setNumberOfIterations(100,100);
    mgsys.solveVcycle();
    //-----------------

    // Make the permanent s_k vector and swap it with the solution
    GenMatrix<UCBspl_real> s_k(n1_,n2_);
    s_k.swap(*s_k_temp); // transfer solution to s_k, but solution vector still refers to the s_k_temp object
 
    //s_k.init(b_);

    GenMatrix<UCBspl_real> p_k;
    p_k.init(s_k);

    // CG p_k.init(r_k);

    // Step 3: rho_0 = r_0^T * r_0
    double rho_k = scalarProduct(s_k,r_k);

    int k=0;
  
    // Step 4 (loop, starting with k=0)
    // --------------------------------  
    GenMatrix<UCBspl_real> t_k(n1_,n2_);
  
    // Precond
    GenMatrix<UCBspl_real> w_k(n1_,n2_);

    while (k < noIterations) {

	// Step 4.1a: t_k = A * p_k
	::matVecMult(A_, p_k, t_k);

	// Step 4.1b: w_k = B * t_k OR run MG cycle on A*w_k = t_k
	// -------------------------------------------------------
	// Note that right hand side changes each time (and so does solution vector)
	// ??? test now and just reinitialize all each time?
	// ??? also change in step 2 later

	s_k_temp->fill(0.0);
	//s_k_temp.init(p_k); // Start vector, but: ??????????? this should set the exact solution !!!!

    
	mgsys.qweSetRHS(t_k); // Set the initial right hand side at the finest level. ??? temporary function now
	mgsys.solveVcycle();

	w_k.swap(*s_k_temp); // transfer solution to w_k
    
	// Step 4.2: alpha_k = rho_k/(p_k^T*t_k)
	double alpha_k = rho_k/scalarProduct(p_k,t_k);

	// Step 4.3: x_(k+1) = x_k + alpha_k * p_k (x_k can be overwritten)
	// (x_k is the algorithm is the private solution vector x_)
	::vecPlusScalVec_1(*x_, p_k, alpha_k);

	// Step 4.4a: r_(k+1) = r_k - alpha_k * t_k (r_k can be overwritten)
	::vecPlusScalVec_1(r_k, t_k, -1.0*alpha_k);

	// Step 4.4b: s_(k+1) = s_k - alpha_k * w_k (s_k can be overwritten)
	::vecPlusScalVec_1(s_k, w_k, -1.0*alpha_k);


	// Step 4.5: rho_(k+1) = s_(k+1)^T * r_(k+1)
	double rho_kPlus1 = scalarProduct(s_k,r_k);


	// Step 4.6: p_(k+1) = s_(k+1) + rho_(k+1)/rho_k * p_k (p_k can be overwritten)
	double rho_fac = rho_kPlus1/rho_k;
	::vecPlusScalVec_2(s_k, p_k, rho_fac);

	// Step 4.7:
	k += 1;

	// Update to next iteration
	// x_k, r_k and p_k have been overwritten and are ok
	rho_k = rho_kPlus1;
    }
}


//------------------------------------------------------------------------
void LSsystem::relaxJacobi(int noIterations, double omega) 
{
    if (noIterations == 0)
	return;

    // Jacobi here:
    GenMatrix<UCBspl_real> x_curr(x_->noX(), x_->noY()); // Unknowns in current iteration

    int j;
    for (int n_it = 1; n_it <= noIterations; n_it++) { 
	for (int i = 0; i < numberOfUnknowns(); i++) { // for each row no. i
	    double sumL = 0.0;
      
	    int offset;
	    for (offset = A_.irow(i); (j=A_.jcol(offset)) < i; offset++) {
		// Find row and column in x (PHI_) from row and column in A
		// (The size of x is  n1 x n2  while the size of A is (n1*n2) x (n1*n2)
        
		int k = j%n1_ - 1; // runs fastest (-1,...,n1-2)
		int l = j/n1_ - 1;
        
		sumL += A_(offset) * (*x_)(k,l);
	    }
      
	    double sumR = 0.0;
      
	    int offset_diag = offset;
	    offset++; // step over diagonal
      
	    for ( ;offset < A_.irow(i+1); offset++) {
		j = A_.jcol(offset);
		int k = j%n1_ - 1; // runs fastest (-1,...,n1-2)
		int l = j/n1_ - 1;
		sumR += A_(offset) * (*x_)(k,l);
	    }
      
	    // Diagonal element in x (and b) is now determined by row no. i in A
	    // Lexicographically or 
	    int k = i%n1_ - 1; // runs fastest (-1,...,n1-2)
	    int l = i/n1_ - 1;

	    // Jacobi here:
	    //x_(k,l) = (b_(k,l) - sumL - sumR)/A_(offset_diag);
	    x_curr(k,l) = (b_(k,l) - sumL - sumR)/A_(offset_diag);
	} // for rows


	// Jacobi here:
	// If omega != 1.0, then weighted or damped Jacobi

	if (omega != 1.0) {
	    int idx = x_->noX()-2;
	    int idy = x_->noY()-2;
	    for (int k = -1; k <= idx; k++)
		for (int l = -1; l <= idy; l++)
		    (*x_)(k,l) = (1.0-omega) * (*x_)(k,l) + omega*(x_curr(k,l));
	}
	else
	    x_->swap(x_curr);
    } // iteration
}

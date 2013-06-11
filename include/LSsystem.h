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

#ifndef _LSSYSTEM_H_
#define _LSSYSTEM_H_

#include <MatSparse.h>
#include <GenMatrix.h>
#include <UCBtypedef.h>

#include <map>

#include <boost/shared_ptr.hpp>

/** \brief <b>Least squares approximation to scattered data with B-splines</b>
 * 
 * LSsystem - The class takes a set of scattered data and produces
 * a B-spline surface that approximates the scattered data.
 * The approximation scheme is a least square fit with
 * a \e thin \e plate \e spline smoothing term.
 * A set of iterative equation solvers can be used to produce the set of
 * B-spline coefficients.
 *
 * No stopping criteria for the equation solvers are provided - the user control
 * is only through the given number of iterations. Thus, we call it \em relaxation 
 * in the documentation.  But the relaxation step may be repeated many times.
 *
 * Note that the format of the scattered data and the solution vector is exactly 
 * the same as in the SINTEF MBA (Multilevel B-spline Approximation) Library. 
 * Thus, a good starting vector can be found by running MBA first. This is also
 *  necessary for surfaces with large coefficient grids (and unevenly distributed
 *  data), which otherwise will converge very slow.
 *
 * This class can also be (and should be) considered as a relaxation step for the
 * multigrid methods in class MGsystem. MGsystem contains a number of multigrid 
 * schemes which produces better results than LSsystem alone.
 *
 * \anchor lssystem_example1
 * Examle of use (minimal):
 * \code
   LSsystem lssys(xArr, yArr, zArr, noX, noY, PHI);  // initialize with scattered data and initial solution vector
   lssys.buildEqSystem(); // Build the equation system
   lssys.relaxCG(noIterations); // Run the given number of iterations with conjugate gradients
   \endcode
 * \author Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
 */
class LSsystem {
  friend class MGsystem;
  

  // Equation system A_ x_ = b_
  MatSparse A_;                                  // The system matrix, with smoothing term included
  boost::shared_ptr<GenMatrix<UCBspl_real> > x_; // The unknown B-spline coefficient matrix (corresponds to PHI_ in MBA)
  GenMatrix<UCBspl_real> b_;                     // Right hand side (also as a matrix) built from z-values

  int n1_, n2_;      // size of spline coefficient matrix: x_.noX(), x_.noY() in A_ x_= b_
  int K_,   L_;      // order of spline surface
  
  double lambda_;    // calculated from the l_2 norms of the matrices: 
  double lambdaFac_; // a factor multiplied with lambda (default = 1)
  
  // The scattered data
  boost::shared_ptr<std::vector<double> > U_;
  boost::shared_ptr<std::vector<double> > V_;
  boost::shared_ptr<std::vector<double> > Zvals_;
  boost::shared_ptr<std::map<int, double, std::less<int> > > weights_;


  std::vector<double> domain_;   // defaults is umin, vmin, umax, vmax
    
  // Internal (private) functions
  void   calculateDomain();
  double umin() const {return domain_[0];}
  double vmin() const {return domain_[1];}
  double umax() const {return domain_[2];}
  double vmax() const {return domain_[3];}
  void   buildSparseStructure();
  void   addSmoothingTerm(double lambdaVal); // (= regularization term) (to A_)
  void   setLimits();
  void   calcLambda(); // for smoothing term

  void buildRHS(); // Build right hand side from data (this is also optionally done in buildEqSystem)

  void residual(GenMatrix<UCBspl_real>& r) const;

public:
  
  /** Constructor with (standard) shared pointers to scattered data
   * \param U: x-values
   * \param V: y-values
   * \param Z: z-values
   * \param noX, noY: Number of coefficients in each direction of the spline surface
   * \param PHI: The spline coefficients
   * \param smoothingFac: Smoothing factor; a weight that determines the smoothness of the surface. See LSsystem::setSmoothingFactor.
   */
  LSsystem(boost::shared_ptr<std::vector<double> > U,
           boost::shared_ptr<std::vector<double> > V,
           boost::shared_ptr<std::vector<double> > Zvals,
           int noX, int noY,       // grid size and number of unknowns is noX*noY
           boost::shared_ptr<GenMatrix<UCBspl_real> > PHI, // solution vector
           double smoothingFac = 1.0);

  ~LSsystem(){}
  
  void setWeights(boost::shared_ptr<std::map<int, double, std::less<int> > > weights) {weights_ = weights;}

  /** Builds the system matrix and the right hand side
   *(uses solution vector given by constructor)
   */
  void buildEqSystem(bool buildRHS = true); // A_ and b_ (and smoothing term)
  
  /** Gauss-Seidel relaxation */
  void relaxGaussSeidel(int noIterations);
  
  /** Conjugate Gradient relaxation */
  void relaxCG(int noIterations);

  /** Conjugate Gradient relaxation with multigrid preconditioning */
  void relaxCGPrecond(int noIterations);

  /** Jacobi relaxation */
  void relaxJacobi(int noIterations, double omega=2./3.);

  /** Set the domain over which the surface is to be defined.
    * The default is the xy-range of the scattered data.
    * If used, this must be done before creating the surface.
    *
    * \note This function can only be used to expand the domain beyond the xy-range
    *       of the scattered data. It is the users responsibility to check that
    *       no scattered data falls outside the domain.
    *       (use std::min_element and std::max_element to find range of data
    *       for std::vector)
    *
    *  \see MGsystem::setDomain
    */
  void setDomain(double umin, double vmin, double umax, double vmax);

  /** Get the domain over which the surface is defined */
  void getDomain(double& umin, double& vmin, double& umax, double& vmax);

  /** Number of unknowns in the equation system, i.e. B-spline coefficients. */
  int  numberOfUnknowns() const {return x_->noX() * x_->noY();} // Which is the same as n1_*n2_

  /** Adjust the thin plate spline energy.
   *  By default this value is 1.0. To get a surface that is smoother (but less accurate)
   *  than the default, a value grater than 1.0 should be given and vice versa.
   *  This will typically be done after relaxation has been run and one want a
   *  a surface that is smoother, or less smooth and more accurate, than obtained in
   *  from previous relaxation. Relaxation must then be run afterwards.
   *
   *  \see MGsystem::setSmoothingFactor
   */
  void setSmoothingFactor(double smoothingFac);

  /** Add a point to the data set.
   *  This can be done after relaxation has been run as in the \ref lssystem_example1 "example" above,
   *  but relaxation must be run again afterwards.
   *  A weight of importance different fro 1.0 can be given.
   *
   ** \retval bool \c false if \c (u,v) is outside the domain.
   */
  bool addPoint(double u, double v, double z, double weight=1.0, bool addToRHS=true);
    
  /** Calculate the l2 norm of the current solution (Frobenius) */
  double currentSolutionNorm() const;

  /** The discrete l_2 norm of the residual ||b-Ax||_2 possibly scaled by area of grid cell */
  double residual_l2(bool scaled=false) const;
  /** The max norm of the residual (unscaled) : ||b-Ax||_inf */
  double residual_linf() const;

  static double rowMatVecMult(int row_no, const MatSparse& A, const GenMatrix<UCBspl_real>& x);
};

#endif


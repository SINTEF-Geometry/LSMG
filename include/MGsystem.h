//================================================================
//
// Created: September 2001
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description: Multigrid scheme
//                                                                           
//================================================================
// Copyright (c) 2002 SINTEF Applied Mathematics
//================================================================
#ifndef _MGSYSTEM_H_
#define _MGSYSTEM_H_

#include <LSsystem.h>
#include <UCBsplineSurface.h>

/** \brief <B>MAIN INTEFFACE: Multigrid schemes for solving least squares approximation to scattered data with B-splines </B>
 *
 * MGsystem - The class takes a set of scattered data as input and produces a
 * B-spline surface that approximates the scattered data.
 * The resulting B-spline surface will be smooth and/or accurate depending
 * on the user input. One of the strengths with the methods implemented in this class
 * is, if desired, smooth and "natural" extrapolation of the surface outside the domain of
 * the scattered data or to areas of the domain with no scattered data.
 *
 * The approximation scheme
 * is a \e least \e square fit with a \e thin \e plate \e spline \e smoothing \e term.
 * A set of multigrid schemes are implemented to solve the equation system that arises.
 * Different choices of relaxation is implemented in the class LSsystem.
 * The functions \ref MGsystem::solveFMG "solveFMG", \ref MGsystem::solveAscend "solveAscend"
 * and \ref MGsystem::solveVcycle "solveVcycle" (and \ref MGsystem::relax "relax")
 * produces the spline surface from the scattered data. This is done by solving a large linear
 * equation system, but due the multigrid schemes this is extremely fast compared
 * to non-multigrid schemes. Though, due to the fact that iterative equation solvers are used,
 * the equation system is not solved exact. But, there are mechanisms for adjusting the initial
 * solution; see below.
 * The spline surface can be retrieved as a SplineSurface object with
 * \ref MGsystem::getSplineSurface "getSplineSurface".
 *
 * No stopping criteria for the equation solvers are provided - the user control is only through
 * the given number of iterations.
 *
 * Note that the format of the scattered data and the solution vector is exactly the same
 * as in the SINTEF MBA (Multilevel B-spline Approximation) Library.
 * Thus, if one uses the V-cycle scheme or simple relaxation only, which starts at the finest
 * level, a good starting vector can be found by running MBA first. This is also necessary for
 * surfaces with large coefficient grids, which otherwise will converge very slowly with the
 * V-cycle scheme or simple relaxation only.
 *
 * \anchor mgsystem_example1
 * Examle of use (minimal with default parameters):
 * \code
 *  MGsystem mgsys();
 *  m0=n0=1; h=7; // Spline surface with 2^h + 3 = 131 coefficients in each direction
 *  mgsys.initMultilevel(m0,n0,h,xArr,yArr,zArr,PHI);
 *  mgsys.solveFMG(); // Solve equation system by running the Full Multigrid Scheme
 *  \endcode
 *
 *  <b>Some hints on usage and setting parameters different from the default</b>.<br>
 *  <UL>
 *  <b>Choice of solver:</b>
 *  The most common multigrid solver is solveFMG, but solveAscend is equally good
 *  in many cases and faster when running with default parameters. These two solvers
 *  should be used if there is no good starting vector present (see above).<br>
 *
 *  <b>Smoothing factor and number of iterations:</b>
 *  The influence of the smoothing factor on the produced B-spline surface
 *  is extensively explained in \ref MGsystem::setSmoothingFactor "setSmoothingFactor".
 *  One should note that when increasing the smoothing factor, more iterations are
 *  needed to obtain a correct surface. This can be controlled by the user in different
 *  ways: One may just set more iterations at each level with
 *  \ref MGsystem::setNumberOfIterations "setNumberOfIterations"; e.g., 50 or 100 iterations
 *  instead of 10 which is the default. But, in most cases this will work better: First run
 *  the solver (for example solveAscend) with the default smoothing factor (1.0).
 *  Then call \ref MGsystem::setSmoothingFactor "setSmoothingFactor"
 *  and run \ref MGsystem::relax "relax". The last command will run the given number
 *  of iterations on the final surface with the new smoothing factor. If you are not
 *  satisfied with the result (see Controlling the result below), you may repeat the
 *  last two steps with new parameters.
 *  Alternatively one may replace relax with \ref MGsystem::solveVcycle "solveVcycle".
 *
 *  <b>Controlling the result and visualizing the surface:</b>
 *  One way of controlling the
 *  result is through \ref MGsystem::residual_l2 "residual_l2" and
 *  \ref MGsystem::residual_linf "residual_linf".
 *  One should note though that even if one or both of the residuals are small, this does not
 *  guarantee that the solution of the equation system is accurate (which is due to an
 *  ill-conditioned equation system).<br>
 *  The safest way to control the quality of the produced B-spline surface is through
 *  visual inspection. The surface can be retrieved at any time
 *  with \ref MGsystem::getSplineSurface "getSplineSurface". This gives a
 *  \ref UCBspl::SplineSurface "SplineSurface" object that can be sampled ( \c z=f(x,y) )
 *  or piped further to utility functions that prints the surface in
 *  visualization format, e.g., VRML.
 * 
 *  <b>Extrapolation:</b>
 *  \anchor mgsysstem_extrapolation 
 *  The mathematical formulation of the B-spline surface used here
 *  (least squares approximation of scattered data with smoothing term)
 *  gives overall good extrapolation properties. That is, the surface can be extrapolated
 *  beyond the xy-range of the input data in a "natural" way by using the function
 *  \ref setDomain before creating the surface.
 *  Such "natural" extrapolation of surfaces is important in many
 *  applications, for example in geological modelling where horizons must be
 *  extrapolated to intersect faults.
 *
 *  </UL>
 *
 * \author Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
 * \see MBA
 */

class MGsystem {

    // Common temporaries for both FMG and V-cycle (but not for solveAscend)
    std::vector<double> domain_;    // defailt is umin, vmin, umax, vmax (see also LSsystem)
    std::vector<LSsystem*> systems_;
    GenMatrix<UCBspl_real> r_;
    GenMatrix<UCBspl_real> e_refined_;
  
    // Array of solution vectors (pointers)
    typedef boost::shared_ptr<GenMatrix<UCBspl_real> > shared_gm;
    std::vector<shared_gm> e_arr_;

    int noIter_;
    int noIterCoarsest_;

    double lambdaFac_;

    int m0_, n0_; // given to init multilevel 
    void calculateStartVector() const;

    double umin() const {return domain_[0];}
    double vmin() const {return domain_[1];}
    double umax() const {return domain_[2];}
    double vmax() const {return domain_[3];}

    void reserveResidual(); // assumes that init has been done

    void residual(const LSsystem& lssys, GenMatrix<UCBspl_real>& r) const;
    void residual_restrictionInjection(const LSsystem& lssys, GenMatrix<UCBspl_real>& r) const;
    void residual_restrictionFullWeightingLaplace(const LSsystem& lssys, 
						  GenMatrix<UCBspl_real>& r) const;
    void residual_restrictionSplines(const LSsystem& lssys, GenMatrix<UCBspl_real>& r) const;

    // Can also be used with a start vector (e.g. to adjust smoothness)
    // ???? Is not using private members now
    void relaxAndCorrectVcycle(int noLevelsInVcycle);

    // Used to control allocation of right hand side (RHS) at the finest level.
    // This is used to avoid building the RHS again if a solveVcycle or relax is
    // run afterwards to improve the solution (or adjust smoothness)
    // false = initial status before any solvers have been run
    // true  = the RHS at the finest level has been calculated
    bool rhsFinestLevel_;

public:
    /** Default parameters should be used through the API. */
    MGsystem() {lambdaFac_=1.0; noIter_=10; noIterCoarsest_=50; rhsFinestLevel_=false;}
    ~MGsystem() {cleanup();}
  
    /** Initialize.
     *
     *  \param m0,n0 (>=1): The initial size of the spline space in the hierarchical construction.
     *               If the rectangular domain is a square, m0=n0=1 is recommended.
     *               If the rectangular domain in the y-direction is twice of that
     *               the x-direction, m0=1, n0=2 is recommended.
     *               In general, if the rectangular domain in the y-direction is 
     *               k times the length in the x-direction, m0=1, n0=k is recommended.
     *               
     *  \param h:    Number of levels in the hierarchical construction.
     *               If, e.g., m0=n0=1 and h=8, The resulting spline surface has a
     *               coefficient grid of size 2^h + 3 = 259 in each direction.
     *  \param U:    Array with x-values of scattered data
     *  \param V:    Array with y-values of scattered data
     *  \param Z:    Array with z-values of scattered data
     *  \param x:    The solution vector, i.e., the tensor product coefficients
     *               (formatted as a matrix).<br>
     *               If the solution vector is allocated to the correct size
     *               and \ref MGsystem::solveVcycle "solveVcycle" or \ref MGsystem::relax "relax" 
     *               is used, the solution vector must also be initialized.
     *               (For example to the average of the z-values, or to zero).
     *
     *  \param weights: (Not required). Each point can be given a weight of importance.
     *         The default weight for each point is 1.0. <br>
     *         Example: Assume that points with index 17 and 80 (starts at zero)
     *         should have weights 5.0 and 2.5 respectively.
     *         The weight argument is then initialized thus:
     *  \code
     *  boost::shared_ptr<std::map<int,double> > weights(new std::map<int,double>);
     *  (*weights.get())[17] = 5.0;
     *  (*weights.get())[80] = 2.5;
     *  \endcode
     */

    typedef boost::shared_ptr<std::map<int, double> > shared_map;
    void initMultilevel(int m0, int n0, int h,
			boost::shared_ptr<std::vector<double> > U,
			boost::shared_ptr<std::vector<double> > V,
			boost::shared_ptr<std::vector<double> > Z,
			boost::shared_ptr<GenMatrix<UCBspl_real> > x, // solution vector
			shared_map weights = shared_map(new std::map<int, double>) );
    //boost::shared_ptr<std::map<int, double> > weights = boost::shared_ptr<std::map<int, double> >(new std::map<int, double>) );

  
              
    /** Full MultiGrid scheme (starting from the finest level). */
    void solveFMG();
  
    /** Simple ascend method from the coarsest to the finest grid (nested iteration) */
    void solveAscend();

    /** Simple ascend method from the coarsest to the finest grid (nested iteration)
     *  that will repeat iterating on each level until the error criterion is less than
     *  'err'. */
    //void solveAscend(double err);

    /** One full V-cycle (starting from the finest level). Run MBA first to find a good starting vector. */
    void solveVcycle();

    /** Run the given number iterations at the finest level.
     * This can be done, e.g., to improve the solution from one of the other solvers,
     * or after a new smoothing parameter has been set with
     * \ref MGsystem::setSmoothingFactor "setSmoothingFactor".
     * See also the detailed description of the class in the header.
     * 
     * \param relaxType 1=Gauss Seidel<br> 2=Jacobi<br> 3=Conjugate Gradients (CG) (recommended)<br> 4=Preconditioned Conjugate gradients
     */
    void relax(int numberOfIterations, int relaxType = 3);

    /** Retrieve the spline surface */
    UCBspl::SplineSurface getSplineSurface() const;

    /** The l_2 norm of the residual ||b - Ax||_2 after solving the equation system,
     * possibly scaled with grid cell size if indicated.
     * A low value indicates a more exact solution of the equation system than a large value.
     */
    double residual_l2(bool scaled=false) const;
    /** The l_inf norm of the residual ||b - Ax||_inf after solving the equation system.
     * A low value indicates a more exact solution of the equation system than a large value.
     */
    double residual_linf() const;

    /** Set the domain over which the surface is to be defined.
     * The default is the xy-range of the scattered data.
     * If used, this must be done before creating the surface.
     *
     * \note This function can only be used to expand the domain beyond the xy-range
     *       of the scattered data, i.e., for \e extrapolation 
     *       (see \ref mgsysstem_extrapolation "above").
     *       It is the users responsibility to check that
     *       no scattered data falls outside the domain.
     *       (use std::min_element and std::max_element to find range of data
     *       for std::vector)
     */
    void setDomain(double umin, double vmin, double umax, double vmax);

    /** Get the domain over which the surface is defined */
    void getDomain(double& umin, double& vmin, double& umax, double& vmax);

    /** Control smoothness of the resulting spline surface. A positive
     * value should be given; default is 1.0.
     *
     * This is a useful function for controlling smoothness versus accuracy of the
     * resulting spline surface.
     * More specifically, the given smoothing factor indicates how strong the
     * thin plate spline energy should influence on the result.
     * A high value gives a surface with "low energy" which is smooth, but it
     * does not necessarily approximate the scattered data well.
     * A low value of the smoothing factor gives a surface which is a good
     * approximation to the scattered data in a least square sense. Thus the
     * accuracy will dominate over smoothness such that the surface may be rough
     * and not pleasant looking since it is forced to approximate all the scattered
     * data well (even noise and outliers are approximated well).
     *
     * 
     * Theoretically, when the smoothing factor increases, the surface converges
     * to a least squares plane through the scattered data.
     * 
     * One should experiment with this function on typical data. For terrain modelling
     * one may try a smoothing factor in the range 0.2 - 1000.0, but lower/higher
     * values may also be useful.
     *
     * \note 
     * A value of \e zero or close to zero may cause break down of the linear equation system
     * in many cases when the scattered data are not uniformly distributed in the domain, and/or
     * the tensor product grid is large compared to the number of input data points.
     * 
     * See also detailed description in the header.
     */
    void setSmoothingFactor(double smoothingFac);

    /** Add a point to the data set.
     *  This can be done after a surface has been calculated from the initial
     *  data set as in the \ref mgsystem_example1 "example" above, and if one want to
     *  add more scattered data points and update the B-spline surface accordingly.
     *  One of the solve...-functions or relax() must be run afterwards
     *  to take the new point into account.
     *  More than one point can be added before solve... and/or relax()
     *  is run again.
     *  The given point is also added to the dataset.
     *  A weight of importance different from 1.0 can be given.
     *
     * \retval bool \c false if \c (u,v) is outside the domain.
     */
    bool addPoint(double u, double v, double z, double weight=1.0);

    /** Set weight of importance of the point with the given index (starting from 0).
     *  Default weight is 1.0
     */
    void setWeight(int pointIndex, double weight);
  
    /** Set number of iterations at each level other than the coarsest. Default is 10 */
    void setNumberOfIterations(int noIter) {noIter_=noIter;}

    /** Set number of iterations at coarsest level.  Default is 50 */
    void setNumberOfIterationsCoarsest(int noIter) {noIterCoarsest_ = noIter;}
  
    /** Clean-up array structures and reduce memory usage */
    void cleanup();

    /* ??? Initialize solution vector with \e value */
    void qweInitX(double value) {
	systems_[systems_.size()-1]->x_->fill(value);
    }

    // ??? temporary function for use in LSsystem::CGPrecond
    /* Set right hand side in equation system. */
    void qweSetRHS(const GenMatrix<UCBspl_real>& t_k) {
	systems_[systems_.size()-1]->b_.init(t_k);
	rhsFinestLevel_ = true;
    }

};

#endif

//! \page mainSimplestMG A complete but simple main program using multigrid schemes in class MGsystem
/// \include mainSimplestMG.cpp
#include <MBA.h>
#include <UCButils.h>
#include <boost/shared_ptr.hpp>
#include <algorithm>

#include <MGsystem.h>


int main() {

  // Read scattered data from file
  // Each array is maintained by "standard" boost shared
  // pointers. (See for example www.boost.org)
  // The format is assumed to be:
  // x y z
  // x y z
  // x y z
  // etc.

  typedef std::vector<double> dVec;
  boost::shared_ptr<dVec> x_arr(new std::vector<double>);
  boost::shared_ptr<dVec> y_arr(new std::vector<double>);
  boost::shared_ptr<dVec> z_arr(new std::vector<double>);
  UCBspl::readScatteredData("Data/scat.dat", *x_arr, *y_arr, *z_arr);
    
  MGsystem mgsys;
  int m0,n0,h;
  m0=n0=1; h=7; // Spline surface with 2^h + 3 = 131 coefficients in
		// each direction
  boost::shared_ptr<GenMatrixType> PHI(new GenMatrixType()); // coefficient matrix (maintained by shared pointer)
  mgsys.initMultilevel(m0,n0,h,x_arr,y_arr,z_arr,PHI);

  // Solve equation system:
  // ----------------------
  mgsys.solveFMG();        // Full MultiGrid scheme
  //mgsys.solveAscend();  // alternatively: Simple ascend method

  // In the case of a simple V-cycle, one should prepare a good
  // starting vector for the coefficient matrix.
  // For example, one may use the MBA algorithm which is extremely
  // fast and produces fair results in most cases;
  // see the other example program for class LSsystem. 
  // mgsys.solveVcycle();  // One V-cycle

  std::cout << "The residuals in l2 and linf norm after FMG: "
	    << mgsys.residual_l2() << " and "
	    << mgsys.residual_linf() << std::endl;

  // Retrieve the spline surface and evaluate
  UCBspl::SplineSurface surface = mgsys.getSplineSurface();  
  
  double xmin = surface.umin(); double ymin = surface.vmin();
  double xmax = surface.umax(); double ymax = surface.vmax();
  double x = (xmin + xmax)/2.; double y = (ymin + ymax)/2.; // Middle of the domain over which the surface is defined
  double nx,ny,nz;

  double z = surface.f(x,y);         // Find height of surface in (x,y).
  surface.normalVector(x,y, nx,ny,nz);      // Find normal vector of
					    // surface in (x,y).

  std::cout << "z-value in (" << x << ',' << y << ") = " <<
      z << std::endl;
  std::cout << "Normal in  (" << x << ',' << y << ") = (" <<
      nx << "," << ny << "," << nz << ")" << std::endl;

  // Sample surface and print to VRML file.
  UCBspl::printVRMLgrid("qwe.wrl", surface, 50, 50, true);
  std::cout << "Printing to qwe.wrl" << std::endl;

  // Make the surface above smoother (at the expence of accuracy)
  mgsys.setSmoothingFactor(100.); // Default was 1.0
  mgsys.relax(100,3);             // 100 iterations with Conjugate
				  // Gradients
  UCBspl::printVRMLgrid("qweSmooth.wrl", surface, 50, 50, true);
  std::cout << "Printing to qweSmooth.wrl" << std::endl;


  return 0;
}


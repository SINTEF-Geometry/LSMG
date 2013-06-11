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


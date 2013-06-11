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

//! \page mainSimplestLS A complete but simple main program using simple relaxation in class LSsystem
/// \include mainSimplestLS.cpp
#include <MBA.h>
#include <UCButils.h>
#include <boost/shared_ptr.hpp>
#include <algorithm>

#include <LSsystem.h>


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

  boost::shared_ptr<GenMatrixType> PHI;
  int noX, noY;
  UCBspl::SplineSurface surface;
  
  bool startVectorFromMBA = true;
  if (startVectorFromMBA) {
    // Optional:
    // Make a good start vector with MBA, The SINTEF Multilevel
    // B-spline Approximation library.
    // See separate documentation for the MBA library
    MBA mba(x_arr, y_arr, z_arr); // Initialize with scattered data
    mba.MBAalg(1,1,7);            // Create spline surface with 2^h +
				  // 3 = 131 coefficients in each
				  // direction
    
    surface = mba.getSplineSurface(); // Retrieve the spline surface
    PHI = mba.PHI(); // The spline coefficient matrix
    
    noX = PHI->noX(); // No. of coefficients in each direction
    noY = PHI->noY();
  }
  else {

    PHI.reset(new GenMatrixType);
    noX = 100; noY = 100;
    PHI->resize(noX, noY);
    // Zoro as starting vector. A better choise may be the mean value
    // of z_vals, or even better, use the MBA library as in the scope
    // above
    PHI->fill(0.0);
    double xmin = *std::min_element(x_arr->begin(), x_arr->end());
    double ymin = *std::min_element(y_arr->begin(), y_arr->end());
    double xmax = *std::max_element(x_arr->begin(), x_arr->end());
    double ymax = *std::max_element(y_arr->begin(), y_arr->end());

    surface.init(PHI, xmin, ymin, xmax, ymax);
  }
  
  LSsystem lssys(x_arr, y_arr, z_arr, noX, noY, PHI);
  lssys.buildEqSystem();
  std::cout << "l2 and linf residual before relaxing: "
	    << lssys.residual_l2() << " and "
	    << lssys.residual_linf() << std::endl;
  int no_iterations = 50;

  lssys.relaxCG(no_iterations); // Relax using Conjugate Gradient method
  std::cout << "l2 and linf residual after relaxing: "
	    <<  lssys.residual_l2() << " and "
	    << lssys.residual_linf() << std::endl;

  // The coefficient matrix PHI is still shared with the
  // UCBspl::SplineSurface object, so we can use it for evaluation.
  
  double x = (surface.umin() + surface.umax())/2.; // In the middle of
						   // the domain
  double y = (surface.vmin() + surface.vmax())/2.;
  double nx,ny,nz;

  double z = surface.f(x,y);         // Find height of surface in (x,y).
  surface.normalVector(x,y, nx,ny,nz);      // Find normal vector of
					    // surface in (x,y).

  std::cout << "z-value in (" << x << ',' << y << ") = "
	    << z << std::endl;
  std::cout << "Normal in  (" << x << ',' << y << ") = ("
	    << nx << "," << ny << "," << nz << ")" << std::endl;

  // Sample surface and print to VRML file.
  UCBspl::printVRMLgrid("qwe.wrl", surface, 50, 50, true);
  std::cout << "Printing to qwe.wrl" << std::endl;

  return 0;
}

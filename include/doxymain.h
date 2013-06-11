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

//! \mainpage Least Squares Approximation of Scattered Data with B-splines
//!
//! \image html Utladalen.jpg "A B-spline surface representing an area in Jotunheimen"
//!
//! <p> 
//! <a name="introLSMG"> <h2>Short introduction</h2> </a>
//! This documentation contains a brief reference manual for the
//! SINTEF <b>LSMG</b> library developed at at <a
//! href="http://www.sintef.no/content/page3____342.aspx">SINTEF Applied
//! Mathematics</a>. The LSMG library is a means to make (penalized)
//! least squares approximations to scattered data with B-spline
//! surfaces.</p>
//! <p>The methods have been developed for making high quality smooth
//! surface approximations to scattered data in applications such as
//! geological modelling, GIS and CAD/CAM. Special attention has been
//! payed to i) creating large surfaces fast and ii) "natural"
//! extrapolation of surfaces outside the domain of the scattered data
//! or to areas of the domain with no scattered data. These properties
//! are important especially in geology and GIS.  Although the
//! approximation methods are global, they are fast due to the
//! multigrid schemes that are used for solving linear equation
//! systems.</p>
//! <p>The main interface is through the class <a class="el"
//! href="classMGsystem.html">MGsystem</a> which contains multigrid
//! schemes for solving the approximation problem fast. The class <a
//! class="el" href="classLSsystem.html">LSsystem</a> which contains
//! basic iterative relaxation procedures (equation solvers) can also
//! be used directly or in combination with <a class="el"
//! href="classMGsystem.html">MGsystem</a>.</p>
//! <p>This library uses the same data structures as the
//! <a href="http://www.sintef.no/upload/IKT/9011/geometri/MBA/mba_doc/index.html">SINTEF MBA library</a>
//! (Multilevel B-spline approximation library) which can optionally
//! be used as a "preconditioner" when running the more CPU demanding
//! algorithms in the LSMG library.</p>
//! <p>
//! <a name="gettingstarted"> <h2>Getting started</h2> </a>
//! It's very easy - just look at the examples referred to below, copy
//! the main programs, modify them for your needs, compile and
//! run!</p>
//! <p>
//! <a name="secMGsystem">
//!   <h2>Solving the scattered data approximation problem with multigrid
//!   equation solvers</h2>
//! </a>
//! These are the most powerful methods for calculating a B-spline
//! surface as a least squares approximation to the scattered
//! data. You can choose among different multigrid schemes, decide the
//! number of iterations to be run by the iterative equation solvers,
//! adjust smoothness of the surface etc. You find a small minimal
//! code <a class="el"
//! href="classMGsystem.html#mgsystem_example1">example</a> in class <a
//! class="el" href="classMGsystem.html">MGsystem</a>. Read the
//! documentation for the functions that are called there. Then, copy
//! the complete <a href="mainSimplestMG.html#mainSimplestMG">main
//! program</a> to a file - compile it and run it! The program samples
//! the resulting spline surface and writes the result to a VRML
//! file.</p>
//! <p>
//! <a name="secLSsystem">
//!   <h2>Simple relaxation</h2>
//! </a>
//! This is a simpler method with pure relaxation only, i.e., starting
//! with a (good) starting vector, use one of the iterative equation
//! solvers in class <a class="el"
//! href="classLSsystem.html">LSsystem</a> and run a given number of
//! iterations. (Here, relaxation means the iteration process by an
//! iterative equation solver). You find a small <a class="el"
//! href="classLSsystem.html#lssystem_example1">example</a> in class <a
//! class="el" href="classLSsystem.html">LSsystem</a>. Read the
//! documentation for the functions that are called there. Then, copy
//! the complete <a href="mainSimplestLS.html#mainSimplestLS">main
//! program</a> to a file - compile it and run it! The program samples
//! the resulting spline surface and writes the result to a
//! VRML-file. This method can be used if you already have a fair
//! initial solution and you want to improve it or smooth it by giving
//! a specified smoothing factor. The initial solution may be the
//! result from the SINTEF MBA library (which is extremely fast), or
//! from the methods in class <a class="el"
//! href="classMGsystem.html">MGsystem</a>.</p>
//! <p>You can download a free VRML-viewer from 
//! <a href="http://www.km.kongsberg.com/sim">Kongsberg SIM</a>.</p>
//! <p>If you have got the Visual C++ workspace with all the source
//! code, then you can build the application and run it with the
//! current main program.</p>
//! <p>
//! <a class="anchor" name="download"></a><h2>Download</h2></p>
//! <p>
//! A GPL-version of the library (for Linux/Unix) can be downloaded
//! from <a
//! href="http://www.sintef.no/math_software">http://www.sintef.no/math_software</a>.
//! <p>Please report any problems or comments to <a
//! href="mailto:jan.b.thomassen@sintef.no">
//! jan.b.thomassen@sintef.no</a></p>
//! <p>Øyvind Hjelle, June 2002</p>
//! <p>
//! <em>Last modified 28.11.2007 by <a
//! href="mailto:jan.b.thomassen@sintef.no">Jan Thomassen</a> </em></p>
//!
//! \image html SINTEFlogo.gif
//!
//! \page mbalib The Multilevel B-splines library
//! <a href="http://www.sintef.no/upload/IKT/9011/geometri/MBA/mba_doc/index.html">Multilevel B-splines library</a>

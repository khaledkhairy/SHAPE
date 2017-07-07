# SHAPE (Spherical HArmonics Parameterization Explorer)

![Alt text](https://github.com/khaledkhairy/SHAPE/blob/master/clips/Screen%20Shot%202016-02-16%20at%2010.05.11%20AM.png "SHAPE screenshot")
-----------------------------------------------------------------------------
Copyright (c) 2017, HHMI-Janelia Research Campus All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  
* Neither the name HHMI-Janelia Research Campus nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL HHMI-JANELIA RESEARCH CAMPUS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
------------------------------------------------------------------------------

## Status: 
In production use at Janelia. This is a nascent set of tools that is undergoing large changes and code cleanup. We consider the library suitable for use by our collaborators as well as other research groups. Due to limited staffing, we do not guarantee support for outside groups.
## Summary: 
Generate spherical harmonics-based shapes, or explore existing ones, and their shape properties interactively and accurately.
SHAPE allows monitoring behavior of (and modification of) spherical harmonics parameterization (SPHARM) coefficients interactively.
These are (Fourier) coefficients of a parametric spherical harmonics representation of objects of genus zero.
Genus zero objects have no handles (the topology is that of the sphere).
SHAPE is also used to view spherical harmonics fields mapped onto a base SPHARM shape.
Shapes generated with SHAPE can be exported as triangular surface meshes into OBJ, PLY, STL and VRML formats.
Test shapes are provided for initial experimentation.

For examples of published works based on shape parameterization using Fourier series expansions, please refer to:

-Minimum-energy vesicle and cell shapes calculated using spherical harmonics parameterization. 
Khaled Khairy and Jonathon Howard. Soft Matter, 2011,7, 2138-2143

-Drawing an elephant with 4 parameters. 
Mayer J., K. Khairy and J. Howard, American Journal of Physics, 78(6):648-649 (2010).

-Shapes of Red Blood Cells: Comparison of 3D Confocal Images with the Bilayer-Couple Model. 
Khaled Khairy, JiJinn Foo and Jonathon Howard. Cellular and Molecular Bioengineering, 1:173-181, 2008

-Detection of Deformable Objects in 3D Images using Markov Chain Monte Carlo and Spherical Harmonics. 
Khairy K., E. Reynaud, E. Stelzer, MICCAI New York 2008

-Spherical harmonics-based parametric deconvolution of 3D surface images using bending energy minimization. 
Khairy Khairy and Jonathon Howard. Med Image Anal. 2008 Apr;12(2):217-27. Epub 2007 Oct 30.

-Spherical Harmonics 3D Active Contours for Membrane Bilayer-Bound Surfaces. 
Khaled Khairy, Jacques Pecreaux and Jonathon Howard. MICCAI 2007

----------------------------------------------------------------------------------
For MacOSX Sierra, it is recommended to use the supplied binary in SHAPE_binary_MACOSX.zip
----------------------------------------------------------------------------------
SHAPE was compiled using Qt 5.7and VTK 7.1.1 (although other versions might also work)

Building notes:

SHAPE has been built using the following C/C++ libraries on both MacOSX and Windows 7:

[1] Eigen 3.1.1: http://eigen.tuxfamily.org

[2] VTK 7.1.1: http://www.vtk.org/

[3] Qt 5.7: http://www.qt.io/

[4] shape tools: https://github.com/khaledkhairy/shape_tools

Building has not been tested extensively on any other configuration.

## Acknowledgements:
To Jonathon Howard for supporting development of the first version of SHAPE.

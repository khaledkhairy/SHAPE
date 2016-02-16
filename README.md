# SHAPE (Spherical HArmonics Parameterization Explorer)

![alt tag] (https://github.com/khaledkhairy/SHAPE/blob/master/clips/Screen%20Shot%202016-02-16%20at%2010.05.11%20AM.png)

Generate spherical harmonics-based shapes, or explore existing ones, and their shape properties interactively and accurately.
SHAPE allows monitoring behavior of (and modification of) spherrical harmonics parameterization (SPHARM) coefficients interactively.
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
For MacOSX, it is recommended to use the supplied binary in SHAPE_binary_MACOSX.zip
----------------------------------------------------------------------------------
SHAPE was compiled using Qt 5.3.0 and VTK 6.10

Building notes:

SHAPE has been built using the following C/C++ libraries on both MacOSX and Windows 7:

[1] Eigen 3.1.1: http://eigen.tuxfamily.org

[2] VTK 6.10: http://www.vtk.org/

[3] Qt 5.3.0: http://www.qt.io/

[4] shape tools: https://github.com/khaledkhairy/shape_tools

Building has not been tested extensively on any other configuration.


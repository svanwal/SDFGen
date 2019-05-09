# About this repository
This code accompanies the paper "The signed distance field (SDF): A faster small-body shape model" by Stefaan Van wal.
The current version is to be used for peer review of the initial version of the manuscript.

This repository includes both C++ code (to construct an SDF from a polyhedron mesh) and Matlab code (to sample the resulting model and run the sample applications discussed in the paper).

The C++ code is based on Christopher Batty's SDFGen tool; the original can be found on GitHub at https://github.com/christopherbatty/SDFGen. The license of this original code and be found in the LICENSE.txt file.

Below, instructions for use of this modified version of SDFGen are given. Separate instructions for the Matlab code can be found in the /Matlab/README.txt file.

# Compiling SDFGen
The SDFGen tool can be easily compiled with cmake. Assuming cmake is installed, you can compile the tool by navigating into the /SDFGen folder in a terminal, and running the following commands:

cmake CMakeLists.txt
make

# Using SDFGen
To use the tool to construct an SDF model, place the source .OBJ file containing the source polyhedron shape model in the /SDFGen folder and run:

./bin/SDFGen <filename> <dx> <padding>

in which:
<filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix .obj
<dx> specifies the length of grid cell in the resulting distance field (in the same units as the source polyhedron).
<padding> specifies the number of cells worth of padding between the object bounding box and the boundary of the distance field grid (must be at least 1).

For example, to construct an SDF of the model <SHAPE_SFM_200k_v20180804.obj> with grid size 2.5 and padding 3, you would run:

./binSDFGen SHAPE_SFM_200k_v20180804.obj 2.5 3

This will generate a .txt file containing the resulting SDF model. This can be read into Matlab using the code and instructions provided in the /Matlab folder of this repository.

# Using different shape models
An .OBJ file of the 200,000-facet model of asteroid Ryugu is included in this repository. The higher-resolution models can be downloaded from http://www.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/ Before using these models, you should use the open-source MeshLab softare (http://www.meshlab.net/) to prepare the models for use:
1) Remove unreferenced vertices (Filters > Cleaning and Repairing > Remove Unreferenced Vertex). You can double-check that the number of vertices and facets (shown on the bottom of the Meshlab window) satisfy Euler's rule as discussed in the paper; SDFGen will perform an independent check of this criterion and let you know if it is not satisfied.
2) Scale the model into proper SI units. The raw models provided by Watanabe et al. are given in kilometers. To scale them into meters, use MeshLab's scaling function (Filters > Normals, Curvatures, and Orientation > Transform: Scale) and enter the value 1000 in the X Axis, Y Axis, and Z Axis fields.
3) Save the mesh (File > Export Mesh As...). Make sure to select .OBJ as the filetype from the dropdown menu. In the window that pops up, make sure to untick the 'Vert_Color' field but keep the 'Vert_Normal' ticked.
The resulting .OBJ file is ready for use with SDFGen.

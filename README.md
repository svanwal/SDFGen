# SDFGen
A simple commandline utility to generate grid-based signed distance field (level set) from polyhedron shape models. This code was adapted by Stefaan Van wal from the original by Christopher Batty (see https://github.com/christopherbatty/SDFGen). The readme of Batty's original code can be found in the LICENSE.txt file.

To compile SDFGen with cmake, run the following commands in a terminal in the /SDFGen folder:
cmake CMakeLists.txt
make

To generate an SDF model, place the source .OBJ file containing the source polyhedron shape model in the /SDFGen folder and run:

./bin/SDFGen <filename> <dx> <padding>

in which:
<filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix .obj
<dx> specifies the length of grid cell in the resulting distance field (in the same units as the source polyhedron).
<padding> specifies the number of cells worth of padding between the object bounding box and the boundary of the distance field grid. Minimum is 1.

This will generate a .txt file containing the resulting SDF model. This can be read into Matlab using the code and instructions provided in the /Matlab folder of this repository. Files for the 200k model of asteroid Ryugu have been provided. The higher-resolution models can be downloaded from http://www.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/

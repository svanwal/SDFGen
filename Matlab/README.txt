License is in a separate file. This file contains the sample Matlab code that runs the examples
from the paper.

It already has a shape model and corresponding SDF of asteroid Ryugu loaded.
If you want to use higher-resolution models, see the notes at the bottom



%% RUNNING NUMERICAL TESTS


%% RUNNING EXAMPLES
Example 1: Run SDFdemo_Altimeter.m replace the files if you want to use the higher res models.
Example 2: Run SDFdemo_Ejecta.m replace the files if you want to use the higher res models.
Example 3: Run SDFdemo_ImageGenFrame.m

%% LOADING MODELS

If you want to use a higher-resolution one, you can download the shape models from
http://www.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/
You should open them with MeshLab which is free
http://www.meshlab.net/
and remove the unreferenced vertices (Filters>Cleaning>Remove unreferenced vertices)
and then scale them to proper units (raw models are in km) (Filters>Normals, Curvatures, and Orientation>Transform: Scale, Normalize)
by 1000 along the X, Y, and Z axes.
Then do File>Export Mesh as... and save it as an .OBJ file.

To construct the SDF, copy the .OBJ shape model file into the folder of SDFGen code,
build SDFGen following the instructions in its readme, and then run
./bin/SDFGen filename dx padding
where dx is the resolution and padding is the number of cells to pad with (3 is fine).

Then copy the resulting .txt file into your Matlab working folder and use the SDFdemo_Loading.m script
to load it into Matlab. Make sure to adjust the filenames for loading and saving correctly.

It uses the loadawobj2016, see the license file in the subfolder. It is unmodified and as provided on the MathWorks website,
see https://www.mathworks.com/matlabcentral/fileexchange/10223-loadawobj

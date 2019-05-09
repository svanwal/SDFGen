# Loading models
The provided Matlab code comes with a loaded version of the SHAPE_SFM_200k_v20180804 model of asteroid Ryugu and a corresponding 5-meter signed distance field. To use different models, they need to be loaded into their proper Matlab structures using the load_models.m file. Simply replace the four filenames at the top of the script (and makes sure Matlab has access to the files).

# Running numerical tests
Two numerical tests are included:
x) Test1_Speed.m: This performs a sampling speed comparison between the polyhedron and SDF, as discussed in Section 3.2 of the paper.
x) Test2_Shape.m: This quantifies the differences between the polyhedron and SDF, as discussed in Section 3.3 of the paper.
These tests can be performed for different shape models by simply replacing the files to be loaded at the top of the scripts.

# Running sample applications
x) Application1_Altimeter.m: Runs the altimeter instrument example of Section 4.1 of the paper.
x) Application2_Ejecta.m: Runs the particle ejecta simulation example of Section 4.2 of the paper.
x) Application3_Shadows.m: Runs the shadow visualization example of Section 4.3 of the paper.
These applications can be performed for different shape models by simply replacing the files to be loaded at the top of the scripts.

# Third-party software
This codes includes a copy of the loadawobj2016 repository by William Harwin. It is unmodified and was obtained from the MathWorks website from https://www.mathworks.com/matlabcentral/fileexchange/10223-loadawobj The license for this code can be found in the /Matlab/loadawobj2016 folder.

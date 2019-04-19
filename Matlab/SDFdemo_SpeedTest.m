% This function loads the OBJ and SDF models that were read by
% SDFdemo_Loading and compares the two shapes
clc
close all
clear all

%% Loading models
disp('Loading models...');
load SHAPE_SFM_200k_v20180804.mat;
load SDF_5m0_SFM_200k_v20180804;

%% Random sampling of the models
Np_sdf = 5e5; % Number of sample points for the SDF
Np_poly = 200; % Number of sample points for the polyhedron
R = 500; % Reference radius of the asteroid
X_sdf = 2*R*rand(3,Np_sdf) - R;
X_poly = 2*R*rand(3,Np_poly) - R;

%% Sampling the SDF
disp(['Sampling ',num2str(Np_sdf),' points with the SDF...']);
tic
for i=1:Np_sdf
%     disp(['SDF Point ',num2str(i),' of ',num2str(Np_sdf)]);
    [d_sdf,N_sdf] = sample_sdf(X_sdf(:,i),sdf);
end
t_sdf = toc;
disp(['The SDF sampling rate was ',num2str(Np_sdf/t_sdf),' distance evaluations per second']);

%% Sampling the polyhedron
disp(['Sampling ',num2str(Np_poly),' points with the polyhedron...']);
tic
for i=1:Np_poly
%     disp(['Polyhedron Point ',num2str(i),' of ',num2str(Np_poly)]);
    [d_poly,N_poly] = sample_poly(X_poly(:,i),poly);
end
t_poly = toc;
disp(['The polyhedron sampling rate was ',num2str(Np_poly/t_poly),' distance evaluations per second']);
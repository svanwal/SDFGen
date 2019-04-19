% This function loads and saves an OBJ model and its SDF in .txt format
clc
close all
clear all

%% Reading files
disp('Loading SDF...');
sdf = read_sdf('SDF_5m0_SFM_200K_v20180804.txt'); % Loads the SDF file
disp('Loading OBJ...');
raw_poly = loadawobj('SHAPE_SFM_200k_v20180804.obj'); % Loads the OBJ file

%% Preparing the polyhedron structure
poly.nPoints = length(raw_poly.v); % Number of vertices
poly.nFacets = length(raw_poly.f3); % Number of facets
poly.tri = raw_poly.f3'; % Facets
poly.pts = raw_poly.v; % Vertices
poly.pN = raw_poly.vn; % Vertex normals
poly.fN = zeros(3,poly.nFacets); % Facet normals
poly.C = zeros(3,poly.nFacets); % Facet centers
for f=1:poly.nFacets % Computing the outward-pointing facet normals and facet centers
    p1 = poly.tri(f,1); P1 = poly.pts(:,p1);
    p2 = poly.tri(f,2); P2 = poly.pts(:,p2);
    p3 = poly.tri(f,3); P3 = poly.pts(:,p3);
    poly.C(:,f) = (P1 + P2 + P3)/3;
    u1 = P2 - P1;
    u2 = P3 - P1;
    N0 = cross(u1,u2);
    if N0'*poly.pN(:,p1)<0
        poly.fN(:,f) = -N0/norm(N0);
    else
        poly.fN(:,f) = N0/norm(N0);
    end
end

%% Saving
save('SDF_5m0_SFM_200k_v20180804.mat','sdf');
save('SHAPE_SFM_200k_v20180804','poly');

%% Plotting
figure(1) % Plots the loaded OBJ polyhedron model
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:));
    axis equal
    shading flat
    colormap bone
    lighting gouraud
    material dull
    light
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    title('Loaded OBJ model');
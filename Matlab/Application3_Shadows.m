% This simulates the shadowing example from Section 4.3 of the paper.
clc
close all
clear all

%% Loading models
disp('Loading models...');
load SHAPE_SFM_200k_v20180804.mat; % The desired polyhedron shape model
load SDF_5m0_SFM_200k_v20180804; % The desired SDF mod

%% Some inputs
X0 = [830;1000;-100]; % Camera location point [m]
uz = [-0.65;-0.75;0.05]; % View direction
phi_sun = deg2rad(8); % Obliquity
n_lambda = 1;
lambda_sun = 3.0; % Time of day in radians relative to prime meridian
b = zeros(1,poly.nFacets); % Brightness

%% Sun location
R_sun = 2000*[cos(phi_sun)*cos(lambda_sun);
              cos(phi_sun)*sin(lambda_sun);
              sin(phi_sun)];

%% Lighting the faces as seen in the 3D figure
u_sun = R_sun/norm(R_sun);
% Determining if the facets are oriented at the sun
facet_angle = u_sun(1)*poly.fN(1,:) + u_sun(2)*poly.fN(2,:) + u_sun(3)*poly.fN(3,:);
id_lit = find(facet_angle>0);
% Raytracing to the SDF intersection along the facet normals
Xs = poly.C(:,id_lit);
Ns = poly.fN(:,id_lit);
d = sample_sdf_multi(Xs,sdf);
disp('Raytracing to intersection...');
for i=1:20
    Xs = Xs - 0.5*[d.*Ns(1,:);d.*Ns(2,:);d.*Ns(3,:)];
    d = sample_sdf_multi(Xs,sdf);
end
% Stepping outside
tol_d = 1e-4;
Ximp = Xs + 10*tol_d*Ns;
% Raytracing to the sun
d = sample_sdf_multi(Ximp,sdf);
disp('Raytracing to sun...');
for i=1:70
    Ximp = Ximp + [u_sun(1)*d;u_sun(2)*d;u_sun(3)*d];
    [d,~,bool_outside] = sample_sdf_multi(Ximp,sdf);
end
% Determining the color
id_lit2 = id_lit(bool_outside==1);
b(id_lit2) = asin(facet_angle(id_lit2));

%% Plotting
figure(1)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:),b(1,:));
    axis equal off
    shading flat
    colormap gray
    set(gca,'clipping','off');
    axis(350*[-1 1 -1 1 -1 1])
    caxis([0 pi/2]);
    view(-150,15);
    axis vis3d
    set(gcf,'color','k');
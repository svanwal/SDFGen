% This function loads the OBJ and SDF models that were read by
% SDFdemo_Loading and compares the two shapes
clc
close all
clear all

%% Loading models
disp('Loading models...');
load SHAPE_SFM_200k_v20180804.mat;
load SDF_5m0_SFM_200k_v20180804.mat;

%% Sampling the SDF at the polyhedron vertices
X_sdf = poly.pts;
tol_d = 1e-4;
disp('Raytracing polyhedron vertices into SDF...');
for p=1:poly.nPoints
    if mod(p,1000)==0
        disp(['   Progress ',num2str(100*p/poly.nPoints),'%']);
    end
    u = poly.pN(:,p)/norm(poly.pN(:,p));
    d = 999;
    while abs(d)>tol_d
        d = sample_sdf(X_sdf(:,p),sdf);
        X_sdf(:,p) = X_sdf(:,p) - 0.9*d*u;
    end
end

%% Sampling the SDF at the polyhedron facet centers
d_sdf = zeros(1,poly.nFacets);
N_sdf = zeros(3,poly.nFacets);
delta_theta = zeros(1,poly.nFacets);
disp('Sampling SDF at polyhedron facets...');
for f=1:poly.nFacets
    [d_sdf(f),N_sdf(:,f)] = sample_sdf(poly.C(:,f),sdf);
    delta_theta(f) = acos(N_sdf(:,f)'*poly.fN(:,f));
end

%% Plotting
figure(1) % Plots the loaded OBJ polyhedron model
subplot(1,2,1)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:));
    axis equal
    shading flat
    colormap bone
    lighting gouraud
    material dull
    light('position',[-1 -1 1]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    title('OBJ Shape');
    set(gca,'clipping','off');
subplot(1,2,2)
    p1 = trisurf(poly.tri,X_sdf(1,:),X_sdf(2,:),X_sdf(3,:));
    axis equal
    shading flat
    colormap bone
    lighting gouraud
    material dull
    light('position',[-1 -1 1]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    title('Raytraced SDF shape');
    set(gca,'clipping','off');
    
figure(2)
subplot(1,2,1)
    p1 = plot(sort(abs(d_sdf)),(1:1:poly.nFacets)/poly.nFacets,'LineWidth',2);
    grid on
    xlabel('|\Delta d| [m]');
    ylabel('cdf [-]');
    set(gca,'fontsize',14);
subplot(1,2,2)
    p1 = plot(sort(rad2deg(delta_theta)),(1:1:poly.nFacets)/poly.nFacets,'LineWidth',2);
    grid on
    xlabel('\Delta \theta [deg]');
    ylabel('cdf [-]');
    set(gca,'fontsize',14);
    
figure(3)
subplot(1,2,1)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:),abs(d_sdf));
    axis equal
    shading flat
    colormap parula
    lighting gouraud
    material dull
    light('position',[-1 -1 1]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    title('Distance error');
    set(gca,'clipping','off');
    c = colorbar;
    ylabel(c,'|\Delta d| [m]','fontsize',14);
    set(c,'fontsize',14);
    caxis([0 0.5]);
subplot(1,2,2)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:),rad2deg(delta_theta));
    axis equal
    shading flat
    colormap parula
    lighting gouraud
    material dull
    light('position',[-1 -1 1]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    title('Normal vector direction error');
    set(gca,'clipping','off');
    c = colorbar;
    ylabel(c,'\Delta \theta [deg]','fontsize',14);
    set(c,'fontsize',14);
    caxis([0 30]);
% This function runs a sample ejecta simulation using the SDF
clc
close all
clear all

%% Loading models
disp('Loading models...');
load SHAPE_SFM_200k_v20180804.mat;
load SDF_5m0_SFM_200k_v20180804.mat;

%% Inputs
Rref = 432.5; % Reference radius
mu = 30; % Gravitational parameter
P = 7.631; % Rotation period in hr
omg = (2*pi)/(3600*P); % Angular velocity
X0 = [1.5*Rref;0;0];
V0 = [0;0;sqrt(mu/norm(X0))] - cross([0;0;omg],X0);

%% Propagating trajectory
dt = 15;
tmax = 5*3600;
opt = odeset('RelTol',1e-8);
info.mu = mu;
info.Omg = [0;0;omg]; 
[t,Y] = ode45(@eom_particle,[0:dt:tmax],[X0;V0],opt,info);
t = t'; Y = Y';

%% Altimeter raytracing to the center
Xs = Y(1:3,:);
alt = zeros(1,numel(t));
tol_d = 1e-4;
for i=1:numel(t)
    d = norm(Y(1:3,i)) - Rref;
    u = -Xs(:,i)/norm(Xs(:,i));
    while abs(d)>tol_d
        d = sample_sdf(Xs(:,i),sdf);
        Xs(:,i) = Xs(:,i) + 0.9*d*u;
    end
    alt(i) = norm(Xs(:,i) - Y(1:3,i));
end
r = sqrt(Y(1,:).^2 + Y(2,:).^2 + Y(3,:).^2);

%% Plotting
df = 5;
figure(1)
cc = colormap(bone(120));
set(gcf,'position',[21 162 978 690]);
    p01 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:));
    hold all
    p04 = plot3([Y(1,1:df:end)' Xs(1,1:df:end)']',[Y(2,1:df:end)' Xs(2,1:df:end)']',[Y(3,1:df:end)' Xs(3,1:df:end)']','-k');
    p02 = plot3(Y(1,:),Y(2,:),Y(3,:),'LineWidth',2);
    p03 = plot3(Xs(1,:),Xs(2,:),Xs(3,:),'LineWidth',2);
    axis equal off
    shading flat
    colormap(cc(20:end,:));
    lighting gouraud
    material dull
    light('position',[1 0 0]);
    view(125,0);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    set(gca,'clipping','off');
    h = legend([p02,p03],' Orbit',' Ground track');
    set(h,'position',[0.4309 0.6272 0.1564 0.0674]);

figure(2)
set(gcf,'position',[1025 329 692 306]);
    p1 = plot(t/3600,alt,'LineWidth',2);
    hold all
    p2 = plot(t/3600,r-Rref,'--','LineWidth',2);
    grid on
    set(gca,'fontsize',14);
    xlabel('t [hr]');
    ylabel('Radial altitude [m]');
    h = legend([p1,p2],'SDF-relative','Sphere-relative','Location','NorthWest');
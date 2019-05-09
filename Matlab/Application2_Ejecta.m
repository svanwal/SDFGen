% This simulates the ejecta example from Section 4.2 of the paper.
clc
close all
clear all

%% Loading models
disp('Loading models...');
load SHAPE_SFM_200k_v20180804.mat; % The desired polyhedron shape model
load SDF_5m0_SFM_200k_v20180804; % The desired SDF model

%% Some inputs
Nsim = 1000; % Number of trajectory simulations
X0 = [34.4;-246.1;369]; % Ejection point [m]
Vmin = 0.10; % Min ejection velocity [m/s]
Vmax = 0.15; % Max ejection velocity [m/s]
alpha_max = deg2rad(15); % Max randomization angle on the impacts
beta_max = deg2rad(5); % Max randomization angle on initial velocity
Omg = 2*pi/(7.631*3600); % Angular velocity of the asteroid [rad/s]
mu = 30; % Gravitational parameter of the asteroid [m3/s2]
eN_min = 0.20; % Min normal restitution coeff
eN_max = 0.50; % Max normal restitution coeff
eT_min = 0.0; % Min tangential restitution coeff
eT_max = 1.0; % Max tangential restitution coeff
v_stop = 0.005; % Post-impact velocity at which the sims are stopped [m/s]

%% Initial conditions
% Making sure the ejection point is outside the asteroid
tol_d = 1e-6;
d0 = 999;
while abs(d0)>tol_d 
    [d0,N0] = sample_sdf(X0,sdf);
    X0 = X0 - d0*N0;
end
X0 = X0 + tol_d*N0;
% Initial positions
X_init = repmat(X0,1,Nsim);
% Initial velocities
V0 = Vmin + (Vmax - Vmin)*rand(1,Nsim);
V_init = [V0*N0(1);V0*N0(2);V0*N0(3)];
beta = beta_max*rand(1,Nsim);
u0 = 2*rand(3,Nsim) - 1;
un = sqrt(u0(1,:).^2 + u0(2,:).^2 + u0(3,:).^2);
u = [u0(1,:)./un;u0(2,:)./un;u0(3,:)./un]; % Random vector
for i=1:Nsim
    u(:,i) = u(:,i) - u(:,i)'*N0*N0;
    u(:,i) = u(:,i)/norm(u(:,i)); % Random vector perpendicular to N0
    V_init(:,i) = myev2dcm(beta(i),u(:,i))*V_init(:,i); % Rotating the velocity vector
end

%% Running the simulations
% Some parameters
dt = 5;
tMax = 5*3600;
max_coll = 25;
pushback = 1e-6;
% Allocating vars
t_imp = nan(Nsim,max_coll);
t_1 = nan(Nsim,1);
t_f = nan(Nsim,1);
X_imp_pre = nan(3,Nsim,max_coll);
V_imp_pre = nan(3,Nsim,max_coll);
X_imp_post = nan(3,Nsim,max_coll);
V_imp_post = nan(3,Nsim,max_coll);
X_1 = nan(3,Nsim);
X_final = nan(3,Nsim);
N_imp_base = nan(3,Nsim,max_coll);
N_imp_rand = nan(3,Nsim,max_coll);
t = cell(1,Nsim);
Y = cell(1,Nsim);
% Defining some options for integration
opt = odeset('RelTol',1e-8,'Events',@event_particle);
info.mu = mu;
info.Omg = [0;0;Omg];
info.sdf = sdf;
% Randomizing the impacts
eN = eN_min + (eN_max - eN_min)*rand(Nsim,max_coll);
eT = eT_min + (eT_max - eT_min)*rand(Nsim,max_coll);
alpha = alpha_max*rand(Nsim,max_coll);
prv0 = 2*rand(3,Nsim,max_coll) - 1;
prv_N = sqrt(prv0(1,:,:).^2 + prv0(2,:,:).^2 + prv0(3,:,:).^2);
prv = zeros(size(prv0));
prv(1,:,:) = prv0(1,:,:)./prv_N;
prv(2,:,:) = prv0(2,:,:)./prv_N;
prv(3,:,:) = prv0(3,:,:)./prv_N;
% Running simulations
n_coll = ones(1,Nsim);
tic;
for i=1:Nsim
    disp(['Propagating trajectory #',num2str(i),' of ',num2str(Nsim)]);
    t{i} = 0;
    Y{i} = [X_init(:,i);V_init(:,i)];
    while n_coll(i)<max_coll+1
%         disp(['   Arc #',num2str(n_coll(i))]);
        % Propagate an arc
        [ttemp,Ytemp] = ode23(@eom_particle,[t{i}(end):dt:tMax],Y{i}(:,end),opt,info);
        % Store the data
        t{i} = [t{i},ttemp'];
        Y{i} = [Y{i},Ytemp'];
        X_imp_pre(:,i,n_coll(i)) = Y{i}(1:3,end);
        V_imp_pre(:,i,n_coll(i)) = Y{i}(4:6,end);
        t_imp(i,n_coll(i)) = t{i}(end);
        t_f(i) = t{i}(end);
        % Impact geometry
        [d_sdf,N_sdf] = sample_sdf(Y{i}(1:3,end),sdf);
        N_imp_base(:,i,n_coll(i)) = N_sdf;
        N_imp_base(:,i,n_coll(i)) = N_imp_base(:,i,n_coll(i)) - N_imp_base(:,i,n_coll(i))'*N_sdf*N_sdf;
        N_imp_base(:,i,n_coll(i)) = N_imp_base(:,i,n_coll(i))/norm(N_imp_base(:,i,n_coll(i)));
        N_imp_rand(:,i,n_coll(i)) = myev2dcm(alpha(i,n_coll(i)),prv(:,i,n_coll(i)))*N_sdf;
        % Evaluating the impact
        X_imp_post(:,i,n_coll(i)) = X_imp_pre(:,i,n_coll(i)) + pushback*N_sdf;
        vN_pre = (V_imp_pre(:,i,n_coll(i))'*N_imp_rand(:,i,n_coll(i)))*N_imp_rand(:,i,n_coll(i));
        vT_pre = V_imp_pre(:,i,n_coll(i)) - vN_pre;
        vN_post = -eN(i,n_coll(i))*vN_pre;
        vT_post = eT(i,n_coll(i))*vT_pre;
        V_imp_post(:,i,n_coll(i)) = vN_post + vT_post;
        % Storing new state
        Y{i}(1:3,end) = X_imp_post(:,i,n_coll(i));
        Y{i}(4:6,end) = V_imp_post(:,i,n_coll(i));
        X_final(:,i) = Y{i}(1:3,end);
        % Increment
        n_coll(i) = n_coll(i) + 1;
        % Checking for settling
        if norm(V_imp_post(:,i,n_coll(i)-1))<v_stop
            break;
        end
    end
    t_1(i) = t_imp(i,1);
    X_1(:,i) = X_imp_pre(:,i,1);
end
tsim = toc;
disp(['The total simulation time for ',num2str(Nsim),' trajectories was ',num2str(tsim),' s.']);
disp(['This is an average performance of ',num2str(Nsim/tsim),' trajectories per second.']);

%% Mvksdensity on settling positions
disp('Evaluating settling point density...');
dX = poly.pts - X0;
r_pts = sqrt(dX(1,:).^2 + dX(2,:).^2 + dX(3,:).^2);
rmax = 300;
bool_active = r_pts<rmax;
f = mvksdensity(X_final',poly.pts(:,bool_active==1)','BandWidth',10);
f2 = zeros(1,poly.nPoints);
f2(bool_active==1) = f;

%% Plotting
figure(1)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:));
    hold all
    p2 = scatter3(X0(1),X0(2),X0(3),20,'filled','markerfacecolor','r');
    p4 = quiver3(X_init(1,1:min(30,Nsim)),X_init(2,1:min(30,Nsim)),X_init(3,1:min(30,Nsim)), ...
                 250*V_init(1,1:min(30,Nsim)),250*V_init(2,1:min(30,Nsim)),250*V_init(3,1:min(30,Nsim)), ...
                 'LineWidth',1,'AutoScale','off','Color','k');
    p3 = quiver3(X0(1),X0(2),X0(3),5*N0(1),5*N0(2),5*N0(3),'LineWidth',2,'Color','r','AutoScale','off');
    for i=1:min(30,Nsim)
        plot3(Y{i}(1,:),Y{i}(2,:),Y{i}(3,:),'-b');
    end
    p5 = scatter3(X_1(1,1:min(30,Nsim)),X_1(2,1:min(30,Nsim)),X_1(3,1:min(30,Nsim)),10,'filled','markerfacecolor','g');
    p5 = scatter3(X_final(1,1:min(30,Nsim)),X_final(2,1:min(30,Nsim)),X_final(3,1:min(30,Nsim)),10,'filled','markerfacecolor','b');
    axis equal off
    shading flat
    colormap bone
    lighting gouraud
    material dull
    light('position',[-1 0 0]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    set(gca,'clipping','on');
    axis([-283.5821  195.5633 -471.2794    8.3747  133.3176  580.9113]);
    view(-40,-3);
    set(gca,'position',[0 0 1 1.2]);
    title('Sample trajectories');
    
figure(2)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:));
    hold all
    p2 = scatter3(X0(1),X0(2),X0(3),20,'filled','markerfacecolor','r');
    p3 = quiver3(X0(1),X0(2),X0(3),25*N0(1),25*N0(2),25*N0(3),'LineWidth',2,'Color','r','AutoScale','off');
    p5 = scatter3(X_final(1,:),X_final(2,:),X_final(3,:),5,'filled','markerfacecolor','b');
    p5 = scatter3(X_1(1,:),X_1(2,:),X_1(3,:),5,'filled','markerfacecolor','g');
    axis equal off
    shading flat
    colormap bone
    lighting gouraud
    material dull
    light('position',[-1 0 0]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    set(gca,'clipping','on');
    axis([-283.5821  195.5633 -471.2794    8.3747  133.3176  580.9113]);
    view(-40,-3);
    set(gca,'position',[0 0 1 1.2]);
    title('Impact and settling points');

figure(3)
    p1 = trisurf(poly.tri,poly.pts(1,:),poly.pts(2,:),poly.pts(3,:),f2);
    hold all
    p2 = scatter3(X0(1),X0(2),X0(3),20,'filled','markerfacecolor','r');
    p3 = quiver3(X0(1),X0(2),X0(3),25*N0(1),25*N0(2),25*N0(3),'LineWidth',2,'Color','r','AutoScale','off');
    axis equal off
    shading interp
    colormap parula
    lighting gouraud
    material dull
    light('position',[-1 0 0]);
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    set(gca,'fontsize',14);
    set(gca,'clipping','on');
    axis([-283.5821  195.5633 -471.2794    8.3747  133.3176  580.9113]);
    view(-40,-3);
    set(gca,'position',[0 0 1 1.2]);
    caxis([0 7e-6]);
    title('Settling point density');
    
figure(4)
    p1 = histogram(t_f/60,'Normalization','pdf','BinWidth',1);
    grid on
    set(gca,'fontsize',14);
    xlabel('Settling time [min]');
    ylabel('pdf [-]');
    title('Settling time statistics');
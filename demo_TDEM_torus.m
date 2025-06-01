% Demo for computing toroidal density-equalizing map of toroidal surfaces 
% using the proposed TDEM method.
%
% Main program:
% [map,uv] = TDEM(v,f,population,R,r,dt,epsilon,max_iter)
% 
% Input:
% v: nv x 3 vertex coordinates of a sliced toroidal surface
% f: nf x 3 triangulations of a sliced toroidal surface
% population: nf x 1 positive quantity
% R: the major radius of the torus
% r: the minor radius of the torus
% dt: step size
% epsilon: stopping parameter
% max_iter: maximum number of iterations
% 
% Output:
% map: nv x 3 vertex coordinates of the toroidal density-equalizing map
% uv: nv x 2 vertex coordinates of the corresponding planar map
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024-2025, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map

addpath('code');
addpath('data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose one example below


load('example_torus_0.mat'); % Illustration example
% load('example_torus_1.mat'); % Example 1
% load('example_torus_2.mat'); % Example 2
% load('example_torus_3.mat'); % Example 3
% load('example_torus_4.mat'); % Example 4
% load('example_torus_5_1.mat'); % Example 5
% load('example_torus_5_2.mat'); % Example 5 (alternative)
% load('example_torus_R2.mat'); % Example with R = 2
% load('example_torus_R4.mat'); % Example with R = 4
% load('example_torus_R6.mat'); % Example with R = 6
% load('example_torus_R8.mat'); % Example with R = 8
% load('example_torus_R10.mat'); % Example with R = 10

% get the glued mesh
[v_ori,f_ori,parent_id] = glue_mesh(v,f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the original torus and the corresponding intiial planar map

% Extract the boundary
bd = meshboundaries(f);
bdy = bd{1};

plot_mesh_with_density(v,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',10);
view([-90 30])
hold on;
% Plot the boundary vertices
plot3(v(bdy([1:end,1]),1),v(bdy([1:end,1]),2),v(bdy([1:end,1]),3),'--','Color', [48 92 222]/255, 'LineWidth',5);
axis equal
title('Initial surface')

% Inverse toroidal projection
uv0 = inverse_toroidal_projection(v,R,r);
plot_mesh_with_density(uv0,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',10);
hold on;
% Plot the boundary vertices
plot(uv0(bdy([1:end,1]),1),uv0(bdy([1:end,1]),2),'--','Color',[48 92 222]/255, 'LineWidth',5);
title('Initial planar domain')

%% Toroidal density-equalizing map
dt = 0.1;
epsilon = 1e-2;
max_iter = 500;

% Main program
tic;
[map,uv] = TDEM(v,f,population,R,r,dt,epsilon,max_iter);
toc;
map_ori = zeros(length(v_ori),3);
map_ori(parent_id,:) = map; % get the final unsliced toroidal map


% Plot the planar mapping result
plot_mesh_with_density(uv,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',3);
hold on;
% Plot the boundary vertices
plot(uv(bdy([1:end,1]),1),uv(bdy([1:end,1]),2),'--','Color', [48 92 222]/255, 'LineWidth',5);
title('Final TDEM planar map')

% Plot the toroidal mapping result
plot_mesh_with_density(map_ori,f_ori,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',3);
view([-90 30])
hold on;
% Plot the boundary vertices
plot3(map(bdy([1:end,1]),1),map(bdy([1:end,1]),2),map(bdy([1:end,1]),3),'--','Color', [48 92 222]/255, 'LineWidth',5);
axis equal
title('Final TDEM toroidal map')

%% Plot the density histograms to assess the mapping performance

density0 = (population/sum(population))./(face_area(f,v)/sum(face_area(f,v)));
density = (population/sum(population))./(face_area(f,map)/sum(face_area(f,map)));

% Find maximum histogram values for consistent histogram axes
counts1 = histcounts(density0,0.05:0.1:3); 
counts2 = histcounts(density,0.05:0.1:3); 
max_value = max([counts1, counts2]);

% Initial density histogram
figure;
histogram(density0,0.05:0.1:3,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([0 3])
ylim([0 max_value*1.1]); 
xlabel('Normalized density')
ylabel('Number of faces')
set(gca,'FontSize',20);
set(gca,'linewidth',3);
title('Initial density')

% Final density histogram 
figure;
histogram(density,0.05:0.1:3,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
xlim([0 3])
ylim([0 max_value*1.1]); 
xlabel('Normalized density')
ylabel('Number of faces')
set(gca,'FontSize',20);
set(gca,'linewidth',3);
title('Final density')
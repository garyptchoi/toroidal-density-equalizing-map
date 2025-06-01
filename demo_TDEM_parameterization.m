% Demo for computing toroidal density-equalizing parameterization of 
% genus-one surfaces using the proposed TDEM method.
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
addpath('data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load a genus-one surface mesh

load('example_bob.mat'); % Bob from Keenan's repository
% load('example_bob_2.mat'); % Bob with another cut path 
% load('example_bracelet.mat'); % raw-bracelet from Thingi10K 
% load('example_twisted.mat'); % twisted model from 180 Wrapped Tubes
% load('example_rocker_arm.mat'); % Rocker arm from INRIA
% load('example_vertebra.mat'); % Vertebra

% plot the original (unsliced) surface 
plot_mesh(v_ori,f_ori);
view([-20 35])
set(gca,'FontSize',20)
title('Initial surface')
hold on
plot3(v_ori(ee([1:end,1]),1),v_ori(ee([1:end,1]),2),v_ori(ee([1:end,1]),3),'-','Color', 'b', 'LineWidth',5);
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slicing the mesh to get a simply connected open surface

% Slice along a set of prescribed edges ee
[f,v,parent_id] = slice_mesh(f_ori,v_ori,ee);

% Get the boundary vertices after slicing
bd = meshboundaries(f);
bdy_index = bd{1};

% Extract the corner vertices (corresponding to the base point)
% Here the corner vertices may not be correctly ordered
vid = find(parent_id == mode(parent_id));
% Reorder them
vid1 = find(bdy_index == vid(1));
vid2 = find(bdy_index == vid(2));
vid3 = find(bdy_index == vid(3));
vid4 = find(bdy_index == vid(4));
corner = bdy_index(sort([vid1,vid2,vid3,vid4]));

% Rearrange the boundary indices
id = find(bdy_index == corner(1));
bdy_index = bdy_index([id:end,1:id-1]);
id1 = 1;
id2 = find(bdy_index == corner(2));
id3 = find(bdy_index == corner(3));
id4 = find(bdy_index == corner(4));

% Ensure that the boundary edges follow the condition of "Width > height"
length1 = sum(sqrt(sum((v(bdy_index((id1+1):id2),:)-...
    v(bdy_index((id1):(id2-1)),:)).^2,2)));
length2 = sum(sqrt(sum((v(bdy_index((id2+1):id3),:)-...
    v(bdy_index((id2):(id3-1)),:)).^2,2)));
if length1 < length2
    corner = corner([2;3;4;1]);
end

%% Initial parameterization with periodic boundary conditions
uv0 = rectangular_map_periodic2(v,f,corner,2*pi*R,2*pi*r);
map0 = toroidal_projection(uv0,R,r);

map0_ori = zeros(size(v_ori));
map0_ori(parent_id,:) = map0;

% Prescribe some population on the mesh
population = face_area(f,v); % defined on the sliced surface mesh
population_ori = face_area(f_ori,v_ori); % defined on the original surface mesh

density0 = population./face_area(f,map0); 

plot_mesh_with_density(v,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',3);
view([-18 36]) 
title('Initial surface')

plot_mesh_with_density(map0,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',3);
view([40 35]) 
title('Initial toroidal paameterization')

plot_mesh_with_density(uv0,f,density0);
colormap spring; 
set(gca,'FontSize',20);
set(gca,'linewidth',3);
title('Initial planar map')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toroidal density-equalizing map
dt = 1e-1;
epsilon = 1e-2;
max_iter = 500;

tic;
% Main program
[map,uv] = TDEM(map0,f,population,R,r,dt,epsilon,max_iter);
map_ori = zeros(size(v_ori));
map_ori(parent_id,:) = map; % get the final unsliced toroidal map
toc;

% Plot the mapping result
plot_mesh(map_ori,f_ori);
view([15 35])
set(gca,'FontSize',20)
title('Final TDEM parameterization result')

%% Plot the initial and final area distortion
d_area0 = area_distortion(v_ori,f_ori,map0_ori);
d_area = area_distortion(v_ori,f_ori,map_ori);

% Initial and final area distortion
disp(['Initial area distortion = ' num2str(mean(abs(d_area0)))])
disp(['Final area distortion = ' num2str(mean(abs(d_area)))])
disp(['Percentage of improvement  = ' num2str(...
    ((mean(abs(d_area0)))-(mean(abs(d_area))))./mean(abs(d_area0)))])

% Find maximum histogram values for consistent histogram axes
counts = histcounts(d_area,-3.9:0.2:3.9); 

figure;
histogram(d_area0,-3.9:0.2:3.9,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
ylim([0 max(counts)*1.1]); 
xlabel('Area distortion')
ylabel('Number of faces')
set(gca,'FontSize',20);
set(gca,'linewidth',3);
title('Initial area distortion')

figure;
histogram(d_area,-3.9:0.2:3.9,'FaceColor',[201 0 22]/255,'FaceAlpha',0.5); 
ylim([0 max(counts)*1.1]); 
xlabel('Area distortion')
ylabel('Number of faces')
set(gca,'FontSize',20);
set(gca,'linewidth',3);
title('Final area distortion')
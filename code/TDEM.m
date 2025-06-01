function [map,uv] = TDEM(v,f,population,R,r,dt,epsilon,max_iter)

% Computing toroidal density-equalizing maps of genus-one toroidal surfaces
% using the proposed TDEM method.
%
% Input:
% v: nv x 3 vertex coordinates of a sliced toroidal surface
% f: nf x 3 triangulations of a sliced toroidal surface
% population: nf x 1 positive quantity
% R: the major radius of the torus
% r: the minor radius of the torus
% dt: step size (default = 0.1)
% epsilon: stopping parameter (default = 1e-2)
% max_iter: maximum number of iterations (default = 500)
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
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map

if nargin < 6
    dt = 0.1;
end

if nargin < 7
    epsilon = 1e-2;
end

if nargin < 8
    max_iter = 500;
end

% Initial flattening
uv0 = inverse_toroidal_projection(v,R,r);

%% Get the four sides of the planar domain for setting the periodic boundary constraints
bd = meshboundaries(f);
bdy_index = bd{1};

% corner: 4 x 1 vertex indices for the four corners of the planar domain, 
%         with anti-clockwise orientation
%               4 - 3
%               |   |
%               1 - 2

corner = zeros(4,1);

% Rearrange the boundary indices
[~,id] = min(abs(uv0(bdy_index,1)-0)+abs(uv0(bdy_index,2)+pi*r));
bdy_index = bdy_index([id:end,1:id-1]);
id1 = 1;
corner(1) = bdy_index(id1);
[~,id2] = min(abs(uv0(bdy_index,1)-2*pi*R)+abs(uv0(bdy_index,2)+pi*r));
corner(2) = bdy_index(id2);
[~,id3] = min(abs(uv0(bdy_index,1)-2*pi*R)+abs(uv0(bdy_index,2)-pi*r));
corner(3) = bdy_index(id3);
[~,id4] = min(abs(uv0(bdy_index,1)-0)+abs(uv0(bdy_index,2)-pi*r));
corner(4) = bdy_index(id4);

if id2 > id3
    % The boundary orientation is wrong, meaning the input f has wrong orientation
    % We correct the orientation and start over
    warning('The input triangulations are with clockwise orientation!');
    f = fliplr(f);
    bd = meshboundaries(f);
    bdy_index = bd{1};
    id = find(bdy_index == corner(1));
    bdy_index = bdy_index([id:end,1:id-1]);
    id1 = 1;
    id2 = find(bdy_index == corner(2));
    id3 = find(bdy_index == corner(3));
    id4 = find(bdy_index == corner(4));
end

bottom = bdy_index(id1:id2);
right = bdy_index(id2:id3);
top = bdy_index(id3:id4);
left = bdy_index([id4:end,id1]);

%% Toroidal density-equalizing map

% Define the initial density based on the toroidal map
map0 = toroidal_projection(uv0,R,r);
rho_f = population./face_area(f,map0); 
rho_v = f2v_area_periodic2(uv0,f,bottom,top,left,right)*rho_f; 
rho_v(flipud(right)) = rho_v(left);
rho_v(flipud(top)) = rho_v(bottom);

uv = uv0;

var_initial = var(rho_v/mean(rho_v));

n = 0;
rho_v_error = std(rho_v)/mean(rho_v);
disp('Step     std(rho)/mean(rho)');  
disp([num2str(n), '        ',num2str(rho_v_error)]);
    
% Iterative scheme
while rho_v_error >= epsilon
    %%
    % Update rho with the periodic boundary conditions enforced
    L = laplace_beltrami_periodic2(uv,f,bottom,top,left,right);
    A = lumped_mass_matrix_periodic2(uv,f,bottom,top,left,right);
    rho_v_temp = (A+dt*L)\(A*rho_v);    
    rho_v_temp(flipud(right)) = rho_v_temp(left);
    rho_v_temp(flipud(top)) = rho_v_temp(bottom);
    
    if sum(sum(isnan(rho_v_temp)))~=0
        break;
    end

    % Update density gradient with the periodic boundary conditions enforced
    grad_rho_temp_f = compute_gradient(uv,f,rho_v_temp);
    grad_rho_temp_v = f2v_area_periodic2(uv,f,bottom,top,left,right)*grad_rho_temp_f;
    grad_rho_temp_v(flipud(right),:) = grad_rho_temp_v(left,:);
    grad_rho_temp_v(flipud(top),:) = grad_rho_temp_v(bottom,:);
    
    % Update the uv map 
    dmap = -[grad_rho_temp_v(:,1) ./ rho_v_temp, grad_rho_temp_v(:,2) ./ rho_v_temp];
    uv = uv + dmap*dt;

    % Ensure bijectivity of the uv map
    mu = beltrami_coefficient(uv0,f,uv);
    if max(abs(mu))>1
        mu(abs(mu)>1) = mu(abs(mu)>1)./abs(mu(abs(mu)>1))*0.9;
        uv = linear_beltrami_solver(uv0,f,mu,bdy_index,uv(bdy_index,:));
    end
     
    % Obtain the toroidal map
    map = toroidal_projection(uv,R,r); 

    % Re-coupling scheme with the periodicity enforced
    rho_f = population./face_area(f,map);
    rho_v = f2v_area_periodic2(uv,f,bottom,top,left,right)*rho_f;
    rho_v(flipud(right)) = rho_v(left);
    rho_v(flipud(top)) = rho_v(bottom);
    
    n = n + 1;

    rho_v_error_new = std(rho_v)/mean(rho_v);
    disp([num2str(n), '        ',num2str(rho_v_error_new)]);

    if n > max_iter 
        break;
    end

    rho_v_error = rho_v_error_new;


end
var_final = var(rho_v/mean(rho_v));

disp(['Var(rho_initial) = ' num2str(var_initial)]);
disp(['Var(rho_final) = ' num2str(var_final)]);

end













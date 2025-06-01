function map = rectangular_map_periodic2(v,f,corner,a,b)
% Mapping a simply connected open mesh onto a rectangular domain 
% [0,a] x [-b/2,b/2] with:
% - the top and bottom boundaries enforced to be periodic in x 
% - the left and right boundaries enforced to be periodic in y
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% corner: 4 x 1 vertex indices for the four corners of the rectangle, with anti-clockwise orientation
% 4 - 3
% |   |
% 1 - 2
% a: The width of the rectangular domain
% b: The height of the rectangular domain
%
% Output:
% map: nv x 2 vertex coordinates of the rectangular parameterization
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map

%% Initial setup

nv = length(v);

if size(v,1) < size(v,2)
    v = v';
end
if size(f,1) < size(f,2)
    f = f';
end
if size(v,2) == 2
    v = [v, zeros(nv,1)];
end

bd = meshboundaries(f);
bdy_index = bd{1};

% rearrange the boundary indices
id = find(bdy_index == corner(1));
bdy_index = bdy_index([id:end,1:id-1]);

id1 = 1;
id2 = find(bdy_index == corner(2));
id3 = find(bdy_index == corner(3));
id4 = find(bdy_index == corner(4));

if id2 > id3
    % the boundary orientation is wrong, meaning the input f has wrong orientation
    % we correct the orientation and start over
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

%% Step 1: Mapping the input mesh onto the unit disk

bdy_length = sqrt((v(bdy_index,1) - v(bdy_index([2:end,1]),1)).^2 + ...
            (v(bdy_index,2) - v(bdy_index([2:end,1]),2)).^2 + ...
            (v(bdy_index,3) - v(bdy_index([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_length),1);

% arc-length parameterization boundary constraint
for i = 2:length(bdy_length)
    for j = 1:i-1
        partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length);
bdy = exp(theta*1i);

% disk harmonic map
M = cotangent_laplacian(v,f);
[mrow,mcol,mval] = find(M(bdy_index,:));
M = M - sparse(bdy_index(mrow),mcol,mval,nv, nv) + ...
        sparse(bdy_index,bdy_index,ones(length(bdy_index),1),nv, nv);
c = zeros(nv,1); 
c(bdy_index) = bdy;
z = M\c;
disk = [real(z),imag(z)]; 

if sum(sum(isnan(disk))) ~= 0
    % use tutte embedding instead
    disk = tutte_map(v,f,bdy_index,bdy); 
end

%% Step 2: Mapping the unit disk to a rectangle

% compute the generalized Laplacian
mu = beltrami_coefficient(disk,f,v);
Ax = generalized_laplacian(disk,f,mu); 
Ay = Ax;

% set the boundary constraints
bx = zeros(nv,1); 
by = zeros(nv,1); 

bottom = bdy_index(id1:id2);
right = bdy_index(id2:id3);
top = bdy_index(id3:id4);
left = bdy_index([id4:end,id1]);

% top and bottom must have the same x-coordinates
Ax(bottom,:) = 0;
Ax = Ax + sparse([bottom;bottom],[bottom;flipud(top)],[ones(length(bottom),1);-ones(length(top),1)],nv,nv);

% left and right must have the same y-coordinates
Ay(left,:) = 0;
Ay = Ay + sparse([left;left],[left;flipud(right)],[ones(length(left),1);-ones(length(right),1)],nv,nv);

% bottom must have y = 0, top must have y = b
Ay([top;bottom],:) = 0;
Ay([top;bottom],[top;bottom]) = diag(ones(length([top;bottom]),1));
by(top) = b;

% left must have x = 0, right must have x = a
Ax([left;right],:) = 0;
Ax([left;right],[left;right]) = diag(ones(length([left;right]),1));
bx(right) = a;

% solve the generalized Laplace equation
map_x = Ax\bx;
map_y = Ay\by;

map = [map_x, map_y-b/2];

end
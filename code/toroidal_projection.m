function map = toroidal_projection(uv,R,r)
% Obtain the vertex coordinates of the torus from the UV coordinates.
% 
% Input:
% uv: nv x 2 UV coordinates
%     - uv(:,1) should be with range 2*pi*R
%     - uv(:,2) should be with range 2*pi*r
% 
% Output:
% map: nv x 3 vertex coordinates of the torus
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map


u = uv(:,1)/R;
v = uv(:,2)/r; 
% depending on how we define the cuts
% if the cut is along the inner hole, then add the -pi factor

map = [(R+r*cos(v)).*cos(u), (R+r*cos(v)).*sin(u), r*sin(v)];


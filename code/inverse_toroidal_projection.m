function uv = inverse_toroidal_projection(v,R,r)
% Obtain the UV coordinates from the vertex coordinates on the torus.
% 
% Input:
% map: nv x 3 vertex coordinates of the torus
% R: major radius
% r: minor radius
% 
% Output:
% uv: nv x 2 UV coordinates
%     - uv(:,1) should be with range 2*pi*R
%     - uv(:,2) should be with range 2*pi*r
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map

X = v(:,1);
Y = v(:,2);
Z = v(:,3);

u = R*asin(Y./sqrt(X.^2+Y.^2)).*(X>=0).*(Y>=0) + ...
    R*(2*pi+asin(Y./sqrt(X.^2+Y.^2))).*(X>=0).*(Y<0) + ...
    R*(pi-asin(Y./sqrt(X.^2+Y.^2))).*(X<0);
    
v = r*asin(Z/r).*(X.^2+Y.^2>=R^2) + ...
    r*(pi-asin(Z/r)).*(X.^2+Y.^2<R^2).*(Z>0) + ...
    -r*(pi+asin(Z/r)).*(X.^2+Y.^2<R^2).*(Z<=0);

uv = [u,v];

    



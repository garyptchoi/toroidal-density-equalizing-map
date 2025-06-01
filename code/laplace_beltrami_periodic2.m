function L = laplace_beltrami_periodic2(v,f,bottom,top,left,right)
% Compute the cotangent Laplacian with the left-right and top-bottom 
% periodic boundary conditions enforced.
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map
nv = length(v);

f1 = f(:,1); 
f2 = f(:,2); 
f3 = f(:,3);

% edge length
l = [sqrt(sum((v(f2,:) - v(f3,:)).^2,2)),...
    sqrt(sum((v(f3,:) - v(f1,:)).^2,2)),...
    sqrt(sum((v(f1,:) - v(f2,:)).^2,2))];
l1 = l(:,1); 
l2 = l(:,2); 
l3 = l(:,3);

% Heron's formula
s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
% cotangent weight
cot12 = (l1.^2 + l2.^2 - l3.^2)./area/4;
cot23 = (l2.^2 + l3.^2 - l1.^2)./area/4; 
cot31 = (l1.^2 + l3.^2 - l2.^2)./area/4; 

% Handle the botoom-top correspondence
idt = find(ismember(f1,top)|ismember(f2,top)|ismember(f3,top));
f1t = f1(idt);
f2t = f2(idt);
f3t = f3(idt);
cot12t = cot12(idt);
cot23t = cot23(idt);
cot31t = cot31(idt);
for i = 1:length(top)
    f1t(f1t==top(i)) = bottom(length(top)+1-i);
    f2t(f2t==top(i)) = bottom(length(top)+1-i);
    f3t(f3t==top(i)) = bottom(length(top)+1-i);
end

% Handle the left-right correspondence
idr = find(ismember(f1,right)|ismember(f2,right)|ismember(f3,right));
f1r = f1(idr);
f2r = f2(idr);
f3r = f3(idr);
cot12r = cot12(idr);
cot23r = cot23(idr);
cot31r = cot31(idr);
for i = 1:length(right)
    f1r(f1r==right(i)) = left(length(right)+1-i);
    f2r(f2r==right(i)) = left(length(right)+1-i);
    f3r(f3r==right(i)) = left(length(right)+1-i);
end

% construct matrix
II = [f1; f2; f2; f3; f3; f1; f1; f2; f3];
JJ = [f2; f1; f3; f2; f1; f3; f1; f2; f3];
V = [-cot12; -cot12; -cot23; -cot23; -cot31; -cot31; ...
    cot12+cot31; cot12+cot23; cot31+cot23]/2;

IIt = [f1t; f2t; f2t; f3t; f3t; f1t; f1t; f2t; f3t];
JJt = [f2t; f1t; f3t; f2t; f1t; f3t; f1t; f2t; f3t];
Vt = [-cot12t; -cot12t; -cot23t; -cot23t; -cot31t; -cot31t; ...
    cot12t+cot31t; cot12t+cot23t; cot31t+cot23t]/2;

IIr = [f1r; f2r; f2r; f3r; f3r; f1r; f1r; f2r; f3r];
JJr = [f2r; f1r; f3r; f2r; f1r; f3r; f1r; f2r; f3r];
Vr = [-cot12r; -cot12r; -cot23r; -cot23r; -cot31r; -cot31r; ...
    cot12r+cot31r; cot12r+cot23r; cot31r+cot23r]/2;

L = sparse([II;IIt;IIr],[JJ;JJt;JJr],[V;Vt;Vr],nv,nv);

end

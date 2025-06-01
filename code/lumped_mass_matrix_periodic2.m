function A = lumped_mass_matrix_periodic2(v,f,bottom,top,left,right)
% compute the lumped mass matrix for FEM Laplacian, with the left-right and 
% top-bottom periodic boundary conditions enforced.
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
f1 = f(:,1); f2 = f(:,2); f3 = f(:,3);

% edge length
l = [sqrt(sum((v(f2,:) - v(f3,:)).^2,2)), sqrt(sum((v(f3,:) - v(f1,:)).^2,2)), sqrt(sum((v(f1,:) - v(f2,:)).^2,2))];
l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);

% Heron's formula
s = (l1 + l2 + l3)*0.5;
area = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
 
% Handle the bottom-top correspondence
idt = find(ismember(f1,top)|ismember(f2,top)|ismember(f3,top));
f1t = f1(idt);
f2t = f2(idt);
f3t = f3(idt);
areat = area(idt);
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
arear = area(idr);
for i = 1:length(right)
    f1r(f1r==right(i)) = left(length(right)+1-i);
    f2r(f2r==right(i)) = left(length(right)+1-i);
    f3r(f3r==right(i)) = left(length(right)+1-i);
end

% construct matrix
II = [f1; f2; f3];
JJ = [f1; f2; f3];
V = [area; area; area]/3;

IIt = [f1t; f2t; f3t];
JJt = [f1t; f2t; f3t];
Vt = [areat; areat; areat]/3;

IIr = [f1r; f2r; f3r];
JJr = [f1r; f2r; f3r];
Vr = [arear; arear; arear]/3;

A = sparse([II;IIt;IIr],[JJ;JJt;JJr],[V;Vt;Vr],nv,nv);


end

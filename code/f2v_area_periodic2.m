function M = f2v_area_periodic2(v,f,bottom,top,left,right)
% Face to vertex interpolation with area weighting, with the left-right 
% and top-bottom periodic boundary conditions enforced.
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
nf = length(f);

if size(v,2) == 2
    v = [v,zeros(nv,1)];
end

% find area
area = face_area(f,v);
f1 = f(:,1);
f2 = f(:,2);
f3 = f(:,3);

% Handle the botoom-top correspondence
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
JJ = [1:nf, 1:nf, 1:nf]';
V = [area; area; area]/3;

IIt = [f1t; f2t; f3t];
JJt = [idt; idt; idt];
Vt = [areat; areat; areat];

IIr = [f1r; f2r; f3r];
JJr = [idr; idr; idr];
Vr = [arear; arear; arear];

M = sparse([II;IIt;IIr],[JJ;JJt;JJr],[V;Vt;Vr],nv,nf);

% normalize
vertex_area_sum = sum(M,2);
[Mrow,Mcol,Mval] = find(M);
M = sparse(Mrow,Mcol,Mval./vertex_area_sum(Mrow),nv,nf);

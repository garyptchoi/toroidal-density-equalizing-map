function [v_ori,f_ori,parent_id] = glue_mesh(v,f)

% Convert a sliced mesh into a glued (unsliced) mesh.
%
% Input:
% v: nv x 3 vertex coordinates of a sliced mesh
% f: nf x 3 faces
%
% Output:
% v_ori: N x 3 vertex coordintes of the glued mesh
% f_ori: nf x 3 faces
% parent_id: nv x 3 indices of the parent vertex ID for all vertices in the
%            sliced mesh
%
% If you use this code in your work, please cite the following paper:
%    S. Yao and G. P. T. Choi,
%    "Toroidal Density-Equalizing Map for Genus-One Surfaces."
%    Preprint, arXiv:2410.16833, 2024.
%
% Copyright (c) 2024, Shunyu Yao and Gary P. T. Choi
%
% https://github.com/garyptchoi/toroidal-density-equalizing-map

nv = size(v,1);

v_ori = v;
f_ori = f;

% extract the boundary vertices
bd = meshboundaries(f);
bdy_index = sort(bd{1}); 

nb = length(bdy_index);
% glue the mesh along the boundary
for i = nb:-1:1
    id = find(abs(v(bdy_index(1:i-1),1)-v(bdy_index(i),1)) + ...
        abs(v(bdy_index(1:i-1),2)-v(bdy_index(i),2)) + ...
        abs(v(bdy_index(1:i-1),3)-v(bdy_index(i),3)) < 1e-6);

    if ~isempty(id)
        % identify this vertex with the vertex with the smallest id
        f_ori(f_ori == bdy_index(i)) = bdy_index(id(1));
           
        % remove the vertex from the vertex list
        v_ori(bdy_index(i),:) = [];
    
        % update the face list change due to the removal
        f_ori(f_ori >= bdy_index(i)) = f_ori(f_ori >= bdy_index(i)) - 1;
    end
end

parent_id = zeros(nv,1);
for i = 1:nv
    [~,id] = min(abs(v_ori(:,1)-v(i,1))+abs(v_ori(:,2)-v(i,2))+abs(v_ori(:,3)-v(i,3)));
    parent_id(i) = id;
end
function plot_mesh(v,f,arg3)
% Plot a mesh.
% 
% Input: 
% v: nv x 3 vertex coordinates
% f: nf x 3 triangulations
% arg3 (optional): nv x 1 quantity defined on vertices

% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
%
% Copyright (c) 2014-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

graph = figure;  
if nargin < 3
    patch('Faces',f,'Vertices',v,'FaceColor',[255, 208, 208]/255,'LineWidth',0.5);
   

else
    patch('Faces',f,'Vertices',v,'FaceColor','flat','FaceVertexCData',arg3,...
        'EdgeColor','none', 'LineWidth',0.5);
    colormap('Copper');
    shading interp;
    set(gcf,'color','w'); 
    
end
axis equal tight off
ax = gca; ax.Clipping = 'off';



assetData = struct('Vertex',v);
setappdata(gca,'AssetData',assetData);

dcm_obj = datacursormode(graph);
set(dcm_obj,'UpdateFcn',@ID_display,'Enable','on')
end

function txt = ID_display(obj,event_obj)

hAxes = get(get(event_obj,'Target'),'Parent');
assetData = getappdata(hAxes,'AssetData');

pos = get(event_obj,'Position');
if length(pos) == 2
    id = vertex_search([pos(1),pos(2)],assetData.Vertex);
    txt = {['(x,y) : (',num2str(pos(1)),',',num2str(pos(2)),')'],...
        ['vertex ID : ',int2str(id)]};
else

    id = vertex_search([pos(1),pos(2),pos(3)],assetData.Vertex);
    txt = {['(x,y,z) : (',num2str(pos(1)),',',num2str(pos(2)),',',num2str(pos(3)),')'],...
        ['vertex ID : ',int2str(id)]};
end
end


function index = vertex_search(XYZ,vertex)

% This function searches vertex indexes for which the coordinates of the
% corresponding vertices are nearest to the input XY / XYZ.
%
% Inputs
% XYZ : (X,Y)/(X,Y,Z) coordinates of ponits in the form of k x 2 / k x 3
% matrices
% vertex : n x 3 vertices coordinates
%
% Output :
%	index : k x 1 vertex indexes of the points (X,Y)/(X,Y,Z) 
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

if size(XYZ,2) ~= 2 && size(XYZ,2) ~= 3
    error('Input feature points must be a kx2 or kx3 matrix');
end


k = size(XYZ,1);
n = size(vertex,1);
index = zeros(k,1);
v = vertex';

switch size(XYZ,2)
    case 2
		[~,index] = min(sqrt((repmat(v(1,:),k,1)-repmat(XYZ(:,1),1,n)).^2 +...
            (repmat(v(2,:),k,1)-repmat(XYZ(:,2),1,n)).^2),[],2);
    case 3
		[~,index] = min(sqrt((repmat(v(1,:),k,1)-repmat(XYZ(:,1),1,n)).^2 +...
            (repmat(v(2,:),k,1)-repmat(XYZ(:,2),1,n)).^2 +...
            (repmat(v(3,:),k,1)-repmat(XYZ(:,3),1,n)).^2),[],2);
end

if length(unique(index)) ~= length(index)
    warning('Some points are sharing the same vertex index found')
end   

end
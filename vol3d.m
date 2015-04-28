function [model] = vol3d(d_value,dz_value, varargin)
%H = VOL3D Volume render 3-D data. 
% VOL3D uses the orthogonal plane 2-D texture mapping technique for 
% volume rending 3-D data in OpenGL. Use the 'texture' option to fine 
% tune the texture mapping technique. This function is best used with
% fast OpenGL hardware.
%
% H = vol3d('CData',data) Create volume render object from input 
%                         3-D data. Use interp3 on data to increase volume
%                         rendering resolution. Returns a struct 
%                         encapsulating the pseudo-volume rendering object.
%                         XxYxZ array represents scaled colormap indices.
%                         XxYxZx3 array represents truecolor RGB values for
%                         each voxel (along the 4th dimension).
%
% vol3d(...,'Alpha',alpha) XxYxZ array of alpha values for each voxel, in
%                          range [0,1]. Default: data (interpreted as
%                          scaled alphamap indices).
%
% vol3d(...,'Parent',axH)  Specify parent axes. Default: gca.
%
% vol3d(...,'texture','2D')  Only render texture planes parallel to nearest
%                            orthogonal viewing plane. Requires doing
%                            vol3d(h) to refresh if the view is rotated
%                            (i.e. using cameratoolbar).
%
% vol3d(...,'texture','3D')  Default. Render x,y,z texture planes
%                            simultaneously. This avoids the need to
%                            refresh the view but requires faster OpenGL
%                            hardware peformance.
%
% vol3d(H)  Refresh view. Updates rendering of texture planes 
%           to reduce visual aliasing when using the 'texture'='2D'
%           option.
%
% NOTES
% Use vol3dtool (from the original vol3d FEX submission) for editing the
% colormap and alphamap. Adjusting these maps will allow you to explore
% your 3-D volume data at various intensity levels. See documentation on 
% alphamap and colormap for more information.
%
% Use interp3 on input date to increase/decrease resolution of data
%
% Examples:
%
% % Visualizing fluid flow
% v = flow(50);
% h = vol3d('cdata',v,'texture','2D');
% view(3); 
% % Update view since 'texture' = '2D'
% vol3d(h);  
% alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
% 
% % Visualizing MRI data
% load mri.mat
% D = squeeze(D);
% h = vol3d('cdata',D,'texture','3D');
% view(3);  
% axis tight;  daspect([1 1 .4])
% alphamap('rampup');
% alphamap(.06 .* alphamap);
%
% See also alphamap, colormap, opengl, isosurface

% Copyright Joe Conti, 2004
% Improvements by Oliver Woodford, 2008-2009, with permission of the
% copyright holder.

if isstruct(varargin{1})
    model = varargin{1};
    if length(varargin) > 1
       varargin = {varargin{2:end}};
    end
else
    model = localGetDefaultModel;
end


if length(varargin)>1
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'cdata'
            model.cdata = varargin{n+1};
        case 'parent'
            model.parent = varargin{n+1};
        case 'texture'
            model.texture = varargin{n+1};
        case 'alpha'
            model.alpha = varargin{n+1};
      case 'start'
        model.start = varargin{n+1};
        
    end
    
  end
end

if isempty(model.parent)
    model.parent = gca;
end

[model] = local_draw(model,d_value, dz_value);

%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];
model.alpha = [];
model.xdata = [];
model.ydata = [];
model.zdata = [];
model.parent = [];
model.handles = [];
model.texture = '3D';



%------------------------------------------%
function [model,ax] = local_draw(model,d_value, dz_value)

cdata = model.cdata; 
siz = size(cdata);
start = model.start;
 d_value = d_value*1000; % *1000 convert m to mm
 dz_value = dz_value*1000;
% Define [x,y,z]data
if isempty(model.xdata)
      % original code, there is no variables like sizeX, sizeY, sizeZ,
      % model.xdata = [0,siz(2)];
    model.xdata = [d_value*(start(1)-1) (siz(1)+start(1)-1)*d_value];
end
if isempty(model.ydata)
    model.ydata = [d_value*(start(2)-1) (siz(2)+start(2)-1)*d_value];
end
if isempty(model.zdata)
    model.zdata = [dz_value*(start(3)-1) (siz(3)+start(3)-1)*dz_value];
end

try
   delete(model.handles);
catch
end

ax = model.parent;
cam_dir = camtarget(ax) - campos(ax);
[m,ind] = max(abs(cam_dir));

opts = {'Parent',ax,'cdatamapping',[],'alphadatamapping',[],'facecolor','texturemap','edgealpha',0,'facealpha','texturemap','tag','vol3d'};

if ndims(cdata) > 3
    opts{4} = 'direct';
else
    cdata = double(cdata);
    opts{4} = 'scaled';
end

if isempty(model.alpha)
    alpha = cdata;
    if ndims(model.cdata) > 3
        alpha = sqrt(sum(double(alpha).^2, 4));
        alpha = alpha - min(alpha(:));
        alpha = 1 - alpha / max(alpha(:));
        size(aplha)
    end
    
    %opts{6} = 'scaled'; % this one produces transparent margin regions
    opts{6} = 'direct';
else
    alpha = model.alpha;
    if ~isequal(siz(1:3), size(alpha))
        error('Incorrect size of alphamatte');
    end
    
    opts{6} = 'none';
end
%max(alpha)
h = findobj(ax,'type','surface','tag','vol3d');
% for n = 1:length(h)
%   try
%      delete(h(n));
%   catch
%   end
% end

is3DTexture = strcmpi(model.texture,'3D');
handle_ind = 1;
%model.xdata
%model.ydata
%model.zdata

% Create z-slice
if(ind==3 || is3DTexture )    
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(2); model.ydata(1), model.ydata(2)];
  %y = [model.ydata(1), model.ydata(2); model.ydata(1), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  %diff = model.zdata(2)-model.zdata(1);
  %delta = diff/size(cdata,3);
  delta = dz_value;
  %for n = 1:size(cdata,3)
   n = 1;
   cslice = squeeze(cdata(:,:,n,:));
   aslice = double(squeeze(alpha(:,:,n)));
   %h(handle_ind) = surface(x,y,z,cslice,'CData',aslice,opts{:});
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %z = z + delta;
   handle_ind = handle_ind + 1;
  %end
  
   n = size(cdata,3);
   z = z + size(cdata,3)*delta;
   cslice = squeeze(cdata(:,:,n,:));
   aslice = double(squeeze(alpha(:,:,n)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
end

% Create x-slice
if (ind==1 || is3DTexture ) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  %diff = model.xdata(2)-model.xdata(1);
  %delta = diff/size(cdata,2);
  delta = d_value;
  %for n = 1:size(cdata,2)
   n = 1;
   cslice = squeeze(cdata(n,:,:,:));
   aslice = double(squeeze(alpha(n,:,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %h(handle_ind) = mesh(x,y,z,cslice,'alphadata',aslice,opts{:});
   %x = x + delta;
   handle_ind = handle_ind + 1;
  %end
  
  n = size(cdata,1);
  x = x + size(cdata,1)*delta;
  cslice = squeeze(cdata(n,:,:,:));
  aslice = double(squeeze(alpha(n,:,:)));
  h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %h(handle_ind) = mesh(x,y,z,cslice,'alphadata',aslice,opts{:});
   %x = x + delta;
   %handle_ind = handle_ind + 1;
   
end
  
  
% Create y-slice
if (ind==2 || is3DTexture)
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  %diff = model.ydata(2)-model.ydata(1);
  %delta = diff/size(cdata,1);
  delta = d_value;
  %for n = 1:size(cdata,1)
   n = 1;
   cslice = squeeze(cdata(:,n,:,:));
   aslice = double(squeeze(alpha(:,n,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
   %y = y + delta;
   handle_ind = handle_ind + 1;
  %end
  
   n = size(cdata,2);
   y = y + size(cdata,2)*delta;
   cslice = squeeze(cdata(:,n,:,:));
   aslice = double(squeeze(alpha(:,n,:)));
   h(handle_ind) = surface(x,y,z,cslice,'alphadata',aslice,opts{:});
end

model.handles = h;
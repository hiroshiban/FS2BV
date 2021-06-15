function voi=GenerateSphericalVOIfromNATIVEcoords(center_XYZ,radius_mm,voinames,voicolors,save_name,VMRdim)

% Generates spherical VOIs in the NATIVE space.
% function GenerateSphericalVOIfromNATIVEcoords(center_XYZ,radius_mm,voinames,:voicolors,:save_name,:VMRdim)
% (: is optional)
%
% This function generates spheric VOI(s) as the center='center_XYZ' with radius='radius_mm'
%
% [example]
% >> center_XYZ=[-35.3,-80.8,5.0;31.7,-84.7,2.5];
% >> voinames={'KO_left_tal','KO_right_tal'};
% >> save_name='KO_talairach.voi';
% >> voi=GenerateSphericalVOIfromNATIVEcoords(center_XYZ,5,voinames,[],save_name);
%
% [input]
% center_XYZ  : talairach coordinate, num_voi x 3 matrix, [X1,Y1,Z1; X2,Y2,Z2,...]
% radius_mm       : radius of the ROI in mm, [val]
% voinames        : names of VOIs, cell structure, {'ROI_1','ROI_2',...}
% voicolors       : (optional) colors of VOIs, cell structure, {[red1,green1,blue1],[red2,...]}
% save_name       : (optional) if not empty, the generated voi structure will be saved
%                   with this name in the current directory, 'my_vois.voi'
% VMRdim          : (optional) VMR volume dimensions in which the generated VOIs
%                   to be projected. 256 by default (1mm-cubic VMR), use 256 for 0.5mm VMR.
%
% [output]
% voi             : generated VOI structure, BVQX voi data
%
% [dependency]
%
% 1. BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
%
% Created    : "2020-03-30 17:06:18 ban"
% Last Update: "2020-03-30 17:13:31 ban"

% check the input variables
if nargin<3, help(mfilename()); return; end

if nargin<4 || isempty(voicolors)
  voicolors=cell(size(center_XYZ,1),1);
  tcolors=rainbow(size(center_XYZ,1));
  for ii=1:1:size(center_XYZ,1)
    voicolors{ii}=tcolors(ii,:);
  end
end

if nargin<5, save_name=''; end
if ~isempty(save_name) && ~strcmp(save_name(end-3:end),'.voi')
  save_name=strcat(save_name,'.voi');
end

if nargin<6 || isempty(VMRdim), VMRdim=256; end

if size(center_XYZ,2)~=3
  error('center_XYZ should be [npoints x 3] matrix. check input variable');
end

if size(center_XYZ,1)~=length(voinames)
  error('the numbers of center_XYZ and voinames mismatch. check input variable');
end

if ~iscell(voinames), voinames={voinames}; end
if ~iscell(voicolors), voicolors={voicolors}; end

scaling=VMRdim/256;

% display message
fprintf('\n');
if size(center_XYZ,1)==1
  fprintf('Generating spheric VOI (radius: %02dmm) in the NATIVE coordinates\n',radius_mm);
else
  fprintf('Generating spheric VOIs (radius: %02dmm) in the NATIVE coordinates\n',radius_mm);
end

% convert unit of radius
radius_mm=(scaling*radius_mm)^2;

% convert tal to BVQX's voi coords
fprintf('Initializing voi full coordinates...');

% Initially we are presented with a list of voi coordinates in the native space.
center_voi_XYZ=scaling.*center_XYZ;

% generate 1mm full-space VOI coordinate
% here, 'repmat' more than 1 repetition makes the processing slow
% using 'for' loop partly is faster
voi_full_coords=zeros(VMRdim,VMRdim,VMRdim,3);
coordarray1=repmat((1:VMRdim)'-0.5,1,VMRdim); % -0.5 is required to set the center of voxel
coordarray2=repmat((1:VMRdim)-0.5,VMRdim,1);
for idx=1:1:VMRdim
  voi_full_coords(:,:,idx,1)=coordarray1;
  voi_full_coords(idx,:,:,2)=coordarray1;
  voi_full_coords(idx,:,:,3)=coordarray2;
end
fprintf('done.\n');

% initialize voi object
voi=BVQXfile('new:voi');
voi.FileVersion=4;
%voi.SubjectVOINamingConvention='<VOI>_<SUBJ>';
voi.NrOfVOIs=size(center_voi_XYZ,1);
voi.OriginalVMRResolutionX=256/VMRdim;
voi.OriginalVMRResolutionY=256/VMRdim;
voi.OriginalVMRResolutionZ=256/VMRdim;
voi.OriginalVMRFramingCubeDim=VMRdim;

% setting new VOI details
for ii=1:1:size(center_voi_XYZ,1)

  % display message
  fprintf('processing (NATIVE:% 6.1f,% 6.1f,% 6.1f) as % 14s...',...
          center_XYZ(ii,1),center_XYZ(ii,2),center_XYZ(ii,3),voinames{ii});

  % generate spheric coordinates
  coord1=voi_full_coords(:,:,:,1)-center_voi_XYZ(ii,1);
  coord2=voi_full_coords(:,:,:,2)-center_voi_XYZ(ii,2);
  coord3=voi_full_coords(:,:,:,3)-center_voi_XYZ(ii,3);
  distances=coord1.^2+coord2.^2+coord3.^2;
  idx=find(distances<=radius_mm);

  % note that BVQX's coords [x,y,z] = [y,z,x] in MATLAB
  [XX,YY,ZZ]=ind2sub([VMRdim,VMRdim,VMRdim],idx);

  % set the coordinates to voi object
  for vv=1:1:numel(XX)
    voi.VOI(ii).Voxels(vv,:)=[VMRdim/2-XX(vv),VMRdim/2-YY(vv),VMRdim/2-ZZ(vv)];
  end

  % set the additional parameters to voi object
  voi.VOI(ii).Name=voinames{ii};
  voi.VOI(ii).Color=voicolors{ii};
  voi.VOI(ii).NrOfVoxels=numel(XX);

  fprintf('done.\n');

end

% save the VOI
if ~isempty(save_name)
  voi.SaveAs(fullfile(pwd,save_name));
end

% clean up
if ~nargout
  voi.ClearObject(); clear voi;
end

fprintf('completed.\n');

return


%% subfunction
function map = rainbow(m)
%
%   function map = rainbow(m)
%
%   RAINBOW(M) returns an M-by-3 matrix containing an RAINBOW colormap.
%   RAINBOW, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(rainbow)
%
%   See also GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM, JET,
%   COLORMAP, RGBPLOT, HSV2RGB, RGB2HSV.
%
% March 98H.Yamamoto
% Modified by Hiroshi Ban, March 08 2011

if nargin<1, m = size(get(gcf,'colormap'),1); end
% h = (0:m-1)'/max(m,1);
h = (m-1:-1:0)'/max(m,1);
if isempty(h)
  map = [];
else
  map = floor(255*hsv2rgb([h ones(m,2)]));
end

return

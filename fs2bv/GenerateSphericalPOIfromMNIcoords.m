function poi=GenerateSphericalPOIfromMNIcoords(SRF_file,center_mni_XYZ,radius_mm,poinames,output_coordinate,poicolors,save_name)

% Generates spherical POIs in MNI space.
% function GenerateSphericalPOIfromMNIcoords(SRF_file,center_mni_XYZ,radius_mm,poinames,:output_coordinate,:poicolors,:save_name)
% (: is optional)
%
% This function generates spheric POI(s) as the center='center_mni_XYZ' with radius='radius_mm'
%
% [example]
% >> center_mni_XYZ=[-35.3,-80.8,5.0;31.7,-84.7,2.5];
% >> poinames={'KO_left_mni','KO_right_mni'};
% >> save_name='KO_mni.poi';
% >> SRF_file='/HB/3d/zk09_034_LH_inflated.srf';
% >> poi=GenerateSphericalPOIfromMNIcoords(SRF_file,center_mni_XYZ,5,poinames,'TAL',[],save_name);
%
% [input]
% SRF_file        : target SRF file name with relative path
%                   e.g. '/HB/3D/zk09_034/zk09_034_LH_inflated.srf'
%                   the location should be specified as such
%                   the current directory where this function is
%                   called is the origin.
% center_mni_XYZ  : MNI coordinate, num_voi x 3 matrix, [X1,Y1,Z1; X2,Y2,Z2,...]
% radius_mm       : radius of the ROI in mm, [val]
% poinames        : names of POIs, cell structure, {'ROI_1','ROI_2',...}
% output_coordinate : (optional) type of the coordinate of the output ROI, 'MNI' or 'TAL'
%                   'TAL' by default
% poicolors       : (optional) colors of POIs, cell structure, {[red1,green1,blue1],[red2,...]}
% save_name       : (optional) if not empty, the generated poi structure will be saved
%                   with this name in the current directory, 'my_pois.poi'
%
% [output]
% poi             : generated POI structure, BVQX poi data
%
% [dependency]
% BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
%
% Created    : "2012-07-05 19:19:04 ban"
% Last Update: "2020-03-30 17:24:43 ban"

%% check the input variables
if nargin<4, help(mfilename()); return; end

srffile=fullfile(pwd,SRF_file);
if ~exist(srffile,'file')
  error('SRF_file is not found. check input variable.');
end

if nargin<5 || isempty(output_coordinate), output_coordinate='TAL'; end

if nargin<6 || isempty(poicolors)
  poicolors=cell(size(center_mni_XYZ,1),1);
  tcolors=rainbow(size(center_mni_XYZ,1));
  for ii=1:1:size(center_mni_XYZ,1)
    poicolors{ii}=tcolors(ii,:);
  end
end

if ~iscell(poinames), poinames={poinames}; end
if ~iscell(poicolors), poicolors={poicolors}; end

if nargin<7, save_name=''; end
if ~isempty(save_name) && ~strcmp(save_name(end-3:end),'.poi')
  save_name=strcat(save_name,'.poi');
end

if size(center_mni_XYZ,2)~=3
  error('center_mni_XYZ should be [npoints x 3] matrix. check input variable');
end

if size(center_mni_XYZ,1)~=length(poinames)
  error('the numbers of center_mni_XYZ and poinames mismatch. check input variable');
end

%%% load SRF file
% display message
fprintf('Target SRF   : %s\n',fullfile(pwd,SRF_file));
srf_tmp=BVQXfile(srffile);

% check whether reference file is attached to the original srfs
if isempty(strfind(srffile,'RECO')) && isempty(strfind(srffile,'RECOSM')) && ...
   isempty(strfind(srffile,'smoothwm')) && isempty(strfind(srffile,'white'))
  if isempty(srf_tmp.AutoLinkedMTC)
    error('reference srf file is not attached to the SRF_file! check AutoLinkedMTC');
  else
    % load the original srf
    if ~exist(srf_tmp.AutoLinkedMTC,'file')
      error('Linked RECOSM file is not found. check file');
    end
    srffile=srf_tmp.AutoLinkedMTC; % update srffile
    srf=BVQXfile(srffile);
    tgtcoords=srf.VertexCoordinate;
    nmeshes=srf.NrOfVertices;
    srf.ClearObject(); clear srf;
  end
else
  % use the current coordinates
  tgtcoords=srf_tmp.VertexCoordinate;
  nmeshes=srf_tmp.NrOfVertices;
end

% clean up
srf_tmp.ClearObject(); clear srf_tmp;

%% generate VOIs

% display message
fprintf('\n');
if size(center_mni_XYZ,1)==1
  fprintf('Generating spheric POI (radius: %02dmm) from MNI coordinates\n',radius_mm);
else
  fprintf('Generating spheric POIs (radius: %02dmm) from MNI coordinates\n',radius_mm);
end

% convert the coordinate from MNI to TAL.
if strcmpi(output_coordinate,'TAL'), center_mni_XYZ=mni2tal(center_mni_XYZ); end

% convert unit of radius
radius_mm=radius_mm^2;

% convert mni to BVQX's voi coords
fprintf('Initializing poi full coordinates...');

% Initially we are presented with a list of voi coordinates
% in mniairach space. These coordinates range from -128 to 128.
center_poi_XYZ=[128-center_mni_XYZ(:,1),128-center_mni_XYZ(:,2),128-center_mni_XYZ(:,3)];
center_poi_XYZ=center_poi_XYZ(:,[2,3,1]); % this is to omit a strange surface coordinate problem

% initialize poi object
poi=BVQXfile('new:poi');
poi.NrOfPOIs=size(center_poi_XYZ,1);

fprintf('done.\n');

% setting new POI details
jj=0;
for ii=1:1:size(center_poi_XYZ,1)

  % display message
  fprintf('processing TAL(converted from MNI):% 6.1f,% 6.1f,% 6.1f as % 14s...',...
          center_mni_XYZ(ii,1),center_mni_XYZ(ii,2),center_mni_XYZ(ii,3),poinames{ii});

  % get SRF indexes within the current spherical region
  srfcoords=tgtcoords-repmat(center_poi_XYZ(ii,:),size(tgtcoords,1),1);
  srfdist=sum(srfcoords.^2,2);

  idx=find(srfdist<=radius_mm);
  if isempty(idx)
    fprintf('no vertex on the cortex is included in the sphere. skipping...\n');
    continue
  end
  jj=jj+1;

  % set the coordinates
  poi.POI(jj).Vertices=idx;

  % set the additional parameters to voi object
  poi.POI(jj).Name=poinames{ii};
  poi.POI(jj).InfoTextFile=char('');
  poi.POI(jj).Color=poicolors{ii};
  poi.POI(jj).NrOfVertices=numel(idx);
  %poi.POI(jj).LabelVertex=0;

  fprintf('done.\n');

end

% set the other parameters
poi.FromMeshFile=srffile;
poi.NrOfMeshVertices=nmeshes;

% save the POI
if ~isempty(save_name)
  poi.SaveAs(fullfile(pwd,save_name));
end

% clean up
if ~nargout
  poi.ClearObject(); clear poi;
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

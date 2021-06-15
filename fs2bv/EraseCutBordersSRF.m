function EraseCutBordersSRF(srf_as_input,srf_as_output,center_XYZ,radius_mm,srf_space)

% Erases POI borders drawn on the SRF fle (regions to be erased are specified by TAL coords).
% function EraseCutBordersSRF(srf_as_input,srf_as_output,center_XYZ,:radius_mm,:srf_space)
% (: is optional)
%
% This fuction removes POI borders drawn on the SRF file. The removal region(s)
% can be specified by TAL or native coordinates. All the borders within a spherical
% region with (center,radius)=(center_XYZ,radius_mm) will be removed.
%
% [IMPORTANT NOTEs]
% 1. when you set srf_space to 'TAL' or 'MNI', and srf_as_input is not a RECOSM, but an inflated or
%    flat format, the original RECOSM file should be linked in advance. TAL_coords on the surface
%    vertices arond the regions of the borders you want to erase can get by CTRL+left-click on the cortex.
% 2. The srf_as_input should be converted to direct color surface (gray-scale color one)
%    in advance of running this function.
% 3. The deleted borders are painted by white color. If you want to remove these white
%    lines, please load the surface onto BrainVoyager and run
%    Meshes -> Background and Curvature Colors... -> Curvature/Smooth
%
% [example 1]
% >> center_XYZ=[-35.3,-80.8,5.0]; % to delete around KO
% >> SRF_input='/HB/3d/zk09_034_LH_inflated_cut.srf';
% >> SRF_output='/HB/3d/zk09_034_LH_inflated_cut_2.srf';
% >> EraseCutBordersSRF(SRF_input,SRF_output,center_XYZ,5);
%
% [example 2]
% To edit and modify the existing borders.
% >> EraseCutBordersSRF('HB_zk11_052.eccentricity_LH_inflated_borders.srf','HB_zk11_052.eccentricity_LH_inflated_borders.srf',[180,181,175;190,173,166;244,154,150;252,147,155],10);
% >> % Edit and modify the borders manually on BrainVoyager with opening the cortex, 'HB_zk11_052.eccentricity_LH_inflated_borders.srf',
% >> % in which the borders are partially erased by the function above.
% >> % After editing, we have to overwrite the existing *_cut.srf surface with the new borders by the function below.
% >> CopySRFcolors('HB_zk11_052.eccentricity_LH_inflated_borders.srf','HB_zk11_052.polar_LH_inflated.srf','HB_zk11_052.polar_RH_inflated_cut.srf',[0,0,0]);
%
% [input]
% srf_as_input   : target SRF file name with relative path
%                  e.g. '/HB/3D/zk09_034/zk09_034_LH_inflated.srf'
%                  the location should be specified as such
%                  the current directory where this function is
%                  called is the origin.
% srf_as_output  : output (border removed) SRF file name with relative path
%                  e.g. '/HB/3D/zk09_034/zk09_034_LH_inflated.srf'
%                  the location should be specified as such
%                  the current directory where this function is
%                  called is the origin.
% center_XYZ     : center coordinate(s), num_points x 3 matrix, [X1,Y1,Z1; X2,Y2,Z2;...]
%                  the borders within a spheric region with center_XYZ as a center and
%                  radius_mm as a radius is(are) erased.
% radius_mm      : (optional) radius of the region to be included in mm, [val], 5 by default
% srf_space      : (optional) the coordinate of the input SRF file.
%                  one of 'NATIVE', 'ACPC', 'TAL', or 'MNI'. 'TAL' by default.
%
% [output]
% no output variable, the new border-removed SRF file will be saved as srf_as_output
%
% [dependency]
% 1. BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
%
% Created    : "2013-04-09 03:16:21 banh"
% Last Update: "2018-09-04 18:34:13 ban"

% check input variables
if nargin<3, help(mfilename()); return; end
if nargin<4 || isempty(radius_mm), radius_mm=5; end
if nargin<5 || isempty(srf_space), srf_space='TAL'; end

if size(center_XYZ,2)~=3
  error('center_XYZ should be [num_points x 3 (x,y,z)] matrix. check the input variable.');
end

srffile=fullfile(pwd,srf_as_input);
if ~exist(srffile,'file'), error('srf_as_input is not found. check the input variable.'); end

if ~strcmpi(srf_space,'ACPC') && ~strcmpi(srf_space,'TAL') && ...
   ~strcmpi(srf_space,'MNI') && ~strcmpi(srf_space,'NATIVE')
  error('srf_space should be one of ''NATIVE'', ''ACPC'', ''TAL'', or ''MNI''. check the input variable.');
end

% convert unit of radius
radius_mm2=radius_mm^2;

% display message
[inputpath,inputfname,inputext]=fileparts(srf_as_input);
fprintf('Target SRF   : %s%s\n',inputfname,inputext);
inputsrf=BVQXfile(srf_as_input);
inputcoords=inputsrf.VertexCoordinate;

% check whether reference file is attached to the original srf (for inflated or flat)
% and load the original surface
if strcmpi(srf_space,'TAL') || strcmpi(srf_space,'MNI')
  if ( isempty(strfind(inputfname,'RECO')) && isempty(strfind(inputfname,'RECOSM')) )
    if isempty(inputsrf.AutoLinkedMTC)
      error('reference srf file is not attached to the srf_as_input! check the srf.AutoLinkedMTC');
    else
      if ~exist(inputsrf.AutoLinkedMTC,'file')
        error('Linked RECOSM file is not found. check the attached RECOSM file');
      else
        % replace the input srf coordinates with the original surface ones
        tmp_srf=BVQXfile(inputsrf.AutoLinkedMTC);
        inputcoords=tmp_srf.VertexCoordinate;
        tmp_srf.ClearObject(); clear tmp_srf;
      end
    end
  end
end

% Convert coordinates so that they range 0-255 when we set srf_space = 'TAL' or 'MNI'.
% The conversion is required as the XYZ should range from -128 to 128 in these coordinates.
% Also, note that -0.5 is required to set the center coordinate of each voxel (assuming 1mm cubic)
if strcmpi(srf_space,'TAL') || strcmpi(srf_space,'MNI')
  XYZ=[128-center_XYZ(:,1),128-center_XYZ(:,2),128-center_XYZ(:,3)]-0.5;
  XYZ=XYZ(:,[2,3,1]); % this is to omit a strange surface coordinate problem
else
  XYZ=center_XYZ;
end

% deleting borders on the surface
fprintf('\ndeleting borders...\n');
for ii=1:1:size(XYZ,1)

  % display message
  if strcmpi(srf_space,'NATIVE') || strcmpi(srf_space,'ACPC')
    fprintf('[XYZ:% 6.1f,% 6.1f,% 6.1f] with %.2f mm rad --> ',...
            center_XYZ(ii,1),center_XYZ(ii,2),center_XYZ(ii,3),radius_mm);
  elseif strcmpi(srf_space,'TAL')
    fprintf('[TAL:% 6.1f,% 6.1f,% 6.1f] with %.2f mm rad --> ',...
            center_XYZ(ii,1),center_XYZ(ii,2),center_XYZ(ii,3),radius_mm);
  elseif strcmpi(srf_space,'MNI')
    fprintf('[MNI:% 6.1f,% 6.1f,% 6.1f] with %.2f mm rad --> ',...
            center_XYZ(ii,1),center_XYZ(ii,2),center_XYZ(ii,3),radius_mm);
  end

  % get SRF index within the target sphere
  srfcoords=inputcoords-repmat(XYZ(ii,:),size(inputcoords,1),1);
  srfdist=sum(srfcoords.^2,2);
  idx=find(srfdist<=radius_mm2);
  if isempty(idx)
    fprintf('no vertex on the cortex is included in the sphere. skipping...\n');
    continue
  end

  % get border vertices
  bidx=find(inputsrf.VertexColor(idx,2)==0 & ...
            inputsrf.VertexColor(idx,3)==0 & ...
            inputsrf.VertexColor(idx,4)==0);

  % set no border flag on the surface
  inputsrf.VertexColor(idx(bidx),1)=NaN;
  inputsrf.VertexColor(idx(bidx),2:4)=255;

  fprintf('%d vertices were deleted.\n',numel(idx(bidx)));

end % for ii=1:1:size(XYZ,1)
fprintf('completed deleting borders.\n');

% save the results
fprintf('saving the result as: %s...',srf_as_output);
inputsrf.SaveAs(fullfile(pwd,srf_as_output));
fprintf('done.\n');

% clean up
inputsrf.ClearObject(); clear inputsrf;

return

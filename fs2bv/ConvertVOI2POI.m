function poi=ConvertVOI2POI(voifile,srffile,sample_range,save_flg,VOI_to_use)

% Converts BVQX VOI to BVQX POI file.
% function poi=ConvertVOI2POI(voifile,srffile,:sample_range,:save_flg,:VOI_to_use)
% (: is optional)
%
% This function reads BVQX VOI (and the corresponding SRF) file and converts
% its volumetric structure to BVQX POI file format.
%
% [IMPORTANT]
% 1. Both voifile and the corresponding srffile should be generated in the
%    common space (e.g. Native, ACPC, or Talairach (TAL)), otherwise the
%    coordinates would be calculated incorrectly.
% 2. The generated POIs would not be necessarily the same with the original
%    POIs (from which the VOIs are generated) defined by drawing the borders.
%    Please be careful.
%
% [example]
% >> srffile='HB/3d/zk11_052.avg.3d_final_TAL_LH_RECOSM.srf';
% >> voifile='HB/voi_files/all_ROIs_lh.voi';
% >> sample_range=0.5;
% >> save_flg=1;
% >> VOI_to_use={'V1_lh','V2d_lh','V3d_lh'};
% >> poi=ConvertVOI2POI(voifile,srffile,sample_range,save_flg,VOI_to_use);
%
% [input]
% voifile : BVQX VOI file, relative path format.
%           e.g. '../voi_files/all_ROIs_lh.voi'
%           voifile should be specified so that
%           the current directory where this function
%           is called is the origin.
% srffile : BVQX SRF file, a corresponding SRF file
%           which is used to generate POI file.
%           relative path format.
%           e.g. '../3d/zk11_052.inflated_LH.srf'
%           srffile should be specified so that the current
%           directory where this function is called is the origin.
% sample_range : (optional) samping range in mm.
%           the cortical surface vertices within
%           -sample_range <= vertex_index <= sample_range
%           from each center of the TAL voxels are assigned as POI vertices.
%           0.5 by default.
% save_flg : (optional) whether saving smp file, [0|1]
%           if save_flg is set to 1, POI file will be saved
%           as (voifile_directory)/(voifile_name).poi, 0 by default
% VOI_to_use : (optional) a cell structure, if not empty, meshes corresponding to the
%           VOI_to_use will be extracted. if empty, all the VOIs will be extracted.
%           empty by default.
%
% [output]
% poi : BVQX POI data structure
% if save_flg is set to 1, POI file will be saved as
% (voifile_directory)/(voifile_name).poi
%
%
% Created    : "2016-11-14 15:28:18 ban"
% Last Update: "2020-07-13 11:46:04 ban"

% check input variable
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(save_flg), save_flg=0; end
if nargin<4 || isempty(sample_range), sample_range=0.5; end
if nargin<5 || isempty(VOI_to_use), VOI_to_use={}; end

if ~isempty(VOI_to_use) && ~iscell(VOI_to_use), VOI_to_use={VOI_to_use}; end

voifile=fullfile(pwd,voifile);
srffile=fullfile(pwd,srffile);

if ~exist(voifile,'file'), error('voifile not found. check input variable.'); end
if ~exist(srffile,'file'), error('srffile not found. check input variable.'); end

% display message
[voipath,voiname,voiext]=fileparts(voifile);
fprintf('target VOI        : %s%s\n',voiname,voiext);

[srfpath,srfname,srfext]=fileparts(srffile); %#ok
fprintf('corresponding SRF : %s%s\n',srfname,srfext);

% processing VOI
voi=BVQXfile(voifile);

% check VOI to be used
voicounter=0;
voilist=[];
if ~isempty(VOI_to_use)
  for vv1=1:1:length(VOI_to_use)
    for vv2=1:1:voi.NrOfVOIs
      if strcmp(voi.VOI(vv2).Name,VOI_to_use{vv1})
        voicounter=voicounter+1;
        voilist(voicounter)=vv2; %#ok
      end
    end
  end
  if length(VOI_to_use)~=numel(voilist)
    fprintf('!!WARNING!!: some ROI(s) you specified in POI_to_use not found. check input variable.\n');
  end
else
  voilist=1:voi.NrOfVOIs;
end

% processing SRF
srf=BVQXfile(srffile);
if ~isempty(srf.AutoLinkedMTC)
  if ~exist(srf.AutoLinkedMTC,'file')
    % trying to find the reference surface in the same directory with the input srffile (provision for file movement and the unwilling detaching...)
    srfpath=fileparts(srffile);
    [dummy,fname,fext]=fileparts(srf.AutoLinkedMTC); %#ok
    if exist(fullfile(srfpath,[fname,fext]),'file')
      srf=BVQXfile(srf.AutoLinkedMTC);
    else
      error('can not find the reference surface. check the input srffile.');
    end
  else
    srf=BVQXfile(srf.AutoLinkedMTC);
  end
end
fprintf('\n');

% initializing POI data structure
poi=BVQXfile('new:poi');
poi.FromMeshFile=strcat(srfname,srfext);
poi.NrOfMeshVertices=srf.NrOfVertices;
poi.NrOfPOIs=numel(voilist);
for pp=2:1:numel(voilist), poi.POI(pp)=poi.POI(1); end % initialize all Maps

% processing over VOIs
vv=0;
for vvvv=voilist
  vv=vv+1;

  fprintf('processing: %s...',voi.VOI(vvvv).Name);

  poi.POI(vv).Name=voi.VOI(vvvv).Name;
  poi.POI(vv).Color=voi.VOI(vvvv).Color;

  % gather vertices within each of cubic voxel spaces
  %
  % IMPORTANT NOTE:
  % Please be careful in the order of the voxel indices.
  % Specifically, [srf_x, srf_y, srf_z] corresponds to [voi_y, voi_z, voi_x]
  % The axes terminology follows the internal BrainVoyager format.
  % The mapping to Talairach axes is as follows:
  % BV X front -> back = Y in Tal space
  % BV Y top -> bottom = Z in Tal space
  % BV Z left -> right = X in Tal space
  % to covert the Talairach axes ([-128, 128]) in VOI to SRF(VMR) axes ([0, 256]).
  % In Tal space, Right X +; Left X -.

  poi_idx=zeros(poi.NrOfMeshVertices,1);
  if strcmpi(voi.ReferenceSpace,'TAL') || strcmpi(voi.ReferenceSpace,'MNI')
    for vox_idx=1:1:voi.VOI(vvvv).NrOfVoxels
      poi_idx( (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,2)-sample_range) <= srf.VertexCoordinate(:,1) & ...
                 srf.VertexCoordinate(:,1) <= (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,2)+sample_range) & ...
               (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,3)-sample_range) <= srf.VertexCoordinate(:,2) & ...
                 srf.VertexCoordinate(:,2) <= (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,3)+sample_range) & ...
               (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,1)-sample_range) <= srf.VertexCoordinate(:,3) & ...
                 srf.VertexCoordinate(:,3) <= (voi.OriginalVMRFramingCubeDim/2-voi.VOI(vvvv).Voxels(vox_idx,1)+sample_range) )=1;
    end
  else % unknown, native, or ACPC
    for vox_idx=1:1:voi.VOI(vvvv).NrOfVoxels
      poi_idx( (voi.VOI(vvvv).Voxels(vox_idx,2)-sample_range) <= srf.VertexCoordinate(:,1) & ...
                 srf.VertexCoordinate(:,1) <= (voi.VOI(vvvv).Voxels(vox_idx,2)+sample_range) & ...
               (voi.VOI(vvvv).Voxels(vox_idx,3)-sample_range) <= srf.VertexCoordinate(:,2) & ...
                 srf.VertexCoordinate(:,2) <= (voi.VOI(vvvv).Voxels(vox_idx,3)+sample_range) & ...
               (voi.VOI(vvvv).Voxels(vox_idx,1)-sample_range) <= srf.VertexCoordinate(:,3) & ...
                 srf.VertexCoordinate(:,3) <= (voi.VOI(vvvv).Voxels(vox_idx,1)+sample_range) )=1;
    end
  end
  poi_idx=sort(find(poi_idx==1));

  poi.POI(vv).NrOfVertices=numel(poi_idx);
  poi.POI(vv).Vertices=poi_idx;
  if ~isempty(poi.POI(vv).Vertices(1))
    poi.POI(vv).LabelVertex=max(poi.POI(vv).Vertices(1)-1,1);
  else
    poi.POI(vv).LabelVertex='';
  end

  fprintf('done.\n');
end % for vvvv=voilist

% saving the result
if save_flg
  fprintf('saving: %s...',strcat(voiname,'.poi'));
  poi.SaveAs(fullfile(voipath,strcat(voiname,'.poi')));
  fprintf('done.\n');
end

% clean up
srf.ClearObject(); clear srf;
voi.ClearObject(); clear voi;
if ~nargout, poi.ClearObject(); clear poi; end

return

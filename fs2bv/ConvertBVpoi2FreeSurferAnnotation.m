function [vertices,labels,clut]=ConvertBVpoi2FreeSurferAnnotation(POI_dir,prefix_poi)

% Reads BrainVoyager POI files (*.poi) and converts them to FreeSurfer surface annotation files (*.annot).
% function [vertices,labels,clut]=ConvertBVpoi2FreeSurferAnnotation(POI_dir,:prefix_poi)
% (: is optional)
%
% This function reads BrainVoyager POI files (vertex-based ROI files, *.poi) and converts them to
% FreeSurfer annotation files, *.annot (those files are generally stored in ~/fs_subj_dir/label/ after recon-all).
%
% NOTE:
% To create FreeSurfer label files from the generated annotation files, please use the unix-based
% FreeSurfer command, mri_annotation2label. for details, plese see
% https://surfer.nmr.mgh.harvard.edu/fswiki/mri_annotation2label
%
% For FreeSurfer details, please see
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki
% https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
% https://surfer.nmr.mgh.harvard.edu/fswiki/SubcorticalSegmentation
% https://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011
%
% [NOTE on importing ROIs, segmentation, surfaces etc defined by the other software, such as FSL or FreeSurfer into BrainVoyager]
% To import the ROIs defined outside BrainVoyager, I prepared the functions below.
% ConvertNiftiRoi2BVvoi_ProbThres : Converts NII-format ROI probability map to BrainVoayer VOIs
%                                   with thresholding the map values
% ConvertNiftiRoi2BVvoi_Labels    : Converts NII-format ROI probability map to BrainVoayer VOIs
%                                   using the label lookuptable corresponding to the map ID
% ExtractFSLroi            : Extracts specific value(s) from NII based on XML database
% ExtractFSLroiDirect      : Extracts specific value(s) from NII directly for ROI generations
% ConvertFSLroi2BVvoi      : Converts FSL NII ROIs to BrainVoyager VOIs
% ConvertSPMroi2BVvoi      : Converts SPM NII ROIs to BrainVoyager VOIs
% ConvertsAALroi2BVvoi     : Converts SPM AAL antomical tempolate (NII) ROIs to BrainVoyager VOIs
% ConvertFreeSurferAnnotation2BVpoi   : Converts FreeSurfer surface annotations to BrainVoyager POIs
% ConvertFreeSurferParcellation2BVvoi : Converts FreeSurfer MGZ parcellations to BrainVoyager VOIs
% ConvertFreeSurferMGZ2VMR : Converts FreeSurer MGZ T1/ROI files to BrainVoayer VMRs
% ConvertFreeSurferRibbon2BL2VMR : Converts FreeSurfer ribbon.mgz to BrainVoyager *_{LH|RH}_BL2.vmr
% ConvertFreeSurferSurf2SRF: Converts FreeSurer surface files to BrainVoayer SRFs
% ImportFreeSurfer2BrainVoyager : Imports FreeSurfer-processed files into BrainVoyager
% MaskVMRbyFreeSurferSegmentation       : for general masking purposes. Any *.mgz segmentation result can be used as a
%                                         mask (by default, the parameters are tuned to process wm.seg.mgz as a mask)
% MaskVMRbyFreeSurferSegmentation_ribbon: Specific for applying a mask using the white (and gray) matter segmentation
%                                         result in ribbon.mgz.
%                                         Generally, for surface reconstructions, MaskVMRbyFreeSurferSegmentation_ribbon
%                                         gives the better results.
% ConvertBVpoi2FreeSurferAnnotation   : Converts BrainVoyager POIs to FreeSurfer Annotation files. We can further
%                                       generate label or volume ROI files from the generated annotation files using
%                                       FreeSurfer commands.
%
% Some of the VOIs defined in TAL/MNI spaces are already stored in ~/BVQX_hbtools/VOIs.
% To find the VOIs in which specific XYZ coordinates in TAL/MNI space, please use the function below,
% GetAreaNameFromAtlasVOI               : Returns area candidates, in which the input XYZ coordinate(s)
%                                        is(are) belonging to, based on the pre-defined VOI atlases.
%
% [example]
% >> POI_dir='./DC/poi_files';
% >> [vertex,label,clut]=ConvertBVpoi2FreeSurferAnnotation(POI_dir,'');
%
% [input]
% POI_dir    : Target directory that contains POI file(s)
%              e.g. '/HB/hb14_061/poi_files'
%              Target directory should be specified as a relative path format in which
%              the current directory where this function is called is the origin.
% prefix_poi : (optional) prefix of POI file(s) to be processed, e.g. prefix_poi='*_lh'; empty by default.
%
% [output]
% all the variables below are cell structure (e.g. vertices{number_of_input_poi_files})
% vertices   : vector with values running from 0 to size(vertices)-1
% labels     : lookup of annotation values for corresponding vertex index.
% colortable : structure of annotation data. for details, please see read_annotation and write_annotation
% The generated FreeSurfer annotation file is stored in the same directoy with the input POI as strrep(poi_file_name,'.poi','.annot')
%
% [note on how to set the 'prefix_*' variable]
% prefix_* can be set flexibly as below.
% 1. a string: setting an including prefix (string) alone
%    e.g. prefix_*='_TDTS6.0';
%         --> processes files whose names contain '_TDTS6.0'
% 2. a {1 x N} cell string: setting including prefix (string) arrays
%    e.g. prefix_*={'_TDTS6.0','_TSS5.0mm'};
%         --> processes files whose names contain '_TDTS6.0s' or '_TSS5.0mm'.
% 3. a {2 x N} cell string: setting including/excluding prefix (string) arrays
%    e.g. prefix_*={{'_TDTS6.0s','_TSS5.0mm'};{'THP'}};
%         --> processes files whose names contain '_TDTS6.0s'
%             or '_TSS5.0mm' but do not contain 'THP'.
%         prefix_*={'';{'_TDTS6.0s'}};
%         --> processes files whose names do not contain '_TDTS6.0s'.
%         prefix_*={'_TSS5.0mm';''};
%         --> processes files whose names contain '_TSS5.0mm'.
%
% [dependency]
% 1. freesurfer_matlab_tools
% Matlab tools to read/write FreeSurfer files.
% The tools are copied under ~/freesurfer/matlab/ when you install FreeSurfer.
%
% 2. UNIX command line tools such as gzip are required to use this function.
% Especially, if you want to use this function on Windows, please install UNIX command emulator,
% Cygwin or MSYS2, and set an envitonmental path to the tools in advance.
%
%
% Created    : "2018-06-09 14:34:38 ban"
% Last Update: "2018-06-09 18:18:58 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_poi), prefix_poi=''; end

% check the freesurfer directory and files
if ~exist(fullfile(pwd,POI_dir),'dir')
  error('POI_dir not found. check the input variable.');
end

% set files to be processed
poifiles=GetFiles(fullfile(pwd,POI_dir),'*.poi',prefix_poi);
if isempty(poifiles), error('No poi file found. check the input variables.'); end

% set the path to the tools for handling FreeSurfer files on MATLAB
warning off;
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));
warning on;

% display message
fprintf('Target POI directory : %s\n',fullfile(pwd,POI_dir));
fprintf('\n');

%% processing
vertices=cell(length(poifiles),1);
labels=cell(length(poifiles),1);
clut=cell(length(poifiles),1);
for ii=1:1:length(poifiles)

  % display message
  [poipath,poifname,poiext]=fileparts(poifiles{ii});
  fprintf('loading the target POI file onto MATLAB: %s%s...',poifname,poiext);

  % loading surface annotation file
  poi=BVQXfile(poifiles{ii});

  % initializations
  vertices{ii}=[];
  labels{ii}=[];
  clut{ii}.numEntries=poi.NrOfPOIs;
  clut{ii}.struct_names=cell(poi.NrOfPOIs,1);
  clut{ii}.table=zeros(poi.NrOfPOIs,5);

  % conversion from BrainVoayer POIs to FreeSurfer annotations
  for pp=1:1:poi.NrOfPOIs
    vertices{ii}=[vertices{ii};poi.POI(pp).Vertices-1]; % as FreeSurfer vertex ID starts from 0, we have to subtract 1.
    clut{ii}.struct_names{pp}=poi.POI(pp).Name;
    clut{ii}.orig_tab='BV';

    % note on clut.table structure
    % table: n x 5 matrix
    %        Columns 1,2,3 are RGB values for struct color
    %        Column 4 is a flag (usually 0)
    %        Column 5 is the structure ID, calculated from
    %        R + G*2^8 + B*2^16 + flag*2^24
    R=poi.POI(pp).Color(1);
    G=poi.POI(pp).Color(2);
    B=poi.POI(pp).Color(3);
    flag=0;
    clut{ii}.table(pp,:)=[R,G,B,flag,R+G*2^8+B*2^16+flag*2^24];
    labels{ii}=[labels{ii};repmat(clut{ii}.table(pp,5),[poi.POI(pp).NrOfVertices,1])];
  end

  % fill the missing vertices and labels by 'Unknown'
  miss_vert=setdiff(0:poi.NrOfMeshVertices-1,vertices{ii})';
  if ~isempty(miss_vert)
    clut{ii}.numEntries=clut{ii}.numEntries+1;
    clut{ii}.struct_names(2:end+1)=clut{ii}.struct_names;
    clut{ii}.struct_names{1}='Unknown';
    vertices{ii}=[vertices{ii};miss_vert];
    labels{ii}=[labels{ii};repmat(0,[numel(miss_vert),1])];
    clut{ii}.table(2:end+1,:)=clut{ii}.table;
    clut{ii}.table(1,:)=[25,5,25,0,1639705];
  end

  % sorting
  [dummy,idx]=sort(vertices{ii},'ascend');
  vertices{ii}=vertices{ii}(idx);
  labels{ii}=labels{ii}(idx);

  % generating the corresponding ctags (required for MNE etc)
  fid=fopen(fullfile(poipath,[poifname,'.ctab']),'w');
  if fid==-1, error('can not create *.ctab file. check the data structure and the directory permissions'); end
  unknown_flg=0;
  for pp=1:1:length(clut{ii}.struct_names)
    if strcmpi(clut{ii}.struct_names{pp},'Unknown')
      fprintf(fid,'  0  Unknown                           0   0   0    0\n');
      unknown_flg=1;
    else
      fprintf(fid,'% 3d  %- 31s ',pp-unknown_flg,clut{ii}.struct_names{pp});
      fprintf(fid,'% 3s % 3s % 3s% 5d\n',num2str(clut{ii}.table(pp,1)),num2str(clut{ii}.table(pp,2)),num2str(clut{ii}.table(pp,3)),0);
    end
  end
  fclose(fid);

  fprintf('done.\n');

  % saving
  fprintf('saving: %s...',[poifname,poiext,'.annot']);
  write_annotation(fullfile(poipath,[poifname,'.annot']),vertices{ii},labels{ii},clut{ii});
  fprintf('done.\n');
  fprintf('\n');

  % clean up
  poi.ClearObject(); clear poi;

end % for ii=1:1:length(fsfiles)

fprintf('completed.\n');

% remove the path to the tools for handling FreeSurfer files on MATLAB
rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

return

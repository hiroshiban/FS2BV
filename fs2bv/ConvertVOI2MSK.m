function msk=ConvertVOI2MSK(VOI_dir,msk_resolution,prefix_voi,save_flg,VOI_to_use)

% Converts BrainVoyager *.voi file(s) to *.msk file(s).
% function msk=ConvertVOI2MSK(VOI_dir,:msk_resolution,:prefix_voi,:save_flg,:VOI_to_use)
% (: is optional)
%
% This function converts VOI files to cortex-mask files.
% The generated MSK file(s) can be used for masking VTC time series etc.
%
% [example]
% >> msk=ConvertVOI2MSK('../fMRI_preprocessed/JM/hb14_008/voi_files','combined_*',1);
%
% [input]
% VOI_dir    : Target directory that contains VOI file(s)
%              e.g. '/HB/hb14_061/voi_files'
%              Target directory should be specified as a relative path format in which
%              the current directory where this function is called is the origin.
% msk_resolution : resolution of the msk file(s), one of 1,2,3. 2 by default.
% prefix_voi : (optional) prefix of VOI file(s) to be processed, e.g. prefix_msk='combined*'; empty by default.
% save_flg   : (optional) whether saving the generated MSK file(s). if 1, the MSK file(s)
%              wil be saved like ~/(voi_file_directory)/(Name_of_VOI).voi seccessively.
%              0 by default.
% VOI_to_use : (optional) a cell structure, if not empty, time series corresponding to the
%              VOI_to_use will be extracted. if empty, all the VOI time series will be
%              extracted. empty by default.
%
% [output]
% msk        : a cell structure, msk{number_of_voi + 1}. BrainVoyaer MSK object
%              the final +1 is for all-VOI-combined mask.
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
% 1. BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
%
% Created    : "2014-09-30 16:19:51 ban"
% Last Update: "2018-09-04 18:34:21 ban"

% check input variables
if nargin<1, help(mfilename()); if nargout, msk=[]; end; return; end
if nargin<2 || isempty(msk_resolution), msk_resolution=2; end
if nargin<3 || isempty(prefix_voi), prefix_voi=''; end
if nargin<4 || isempty(save_flg), save_flg=0; end
if nargin<5 || isempty(VOI_to_use), VOI_to_use={}; end

if ~isempty(VOI_to_use) && ~iscell(VOI_to_use), VOI_to_use={VOI_to_use}; end

if isempty(intersect(msk_resolution,[1,2,3]))
  error('msk_resolution should be one of 1,2,3. check input variables.');
end

% set a BrainVoyager constant variable. the parameter below is the default TAL box used in the BrainVoyager.
msk_box=[ 57,  52,  59;  % [XStart, YStart, ZStart]
         231, 172, 197]; % [XEnd  , YEnd  , ZEnd  ]

% get the target VOI file
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if isempty(voifiles), error('no VOI file found. check input variable.'); end

% get voi files from the directories in the VOI_dir
% if multiple files are found, they are temporally combined into one VOI
if length(voifiles)>1
  fprintf('multiple voi files found. combining them into one temporally...\n');
  CombineVOIs(VOI_dir,'tmp_combined.voi',prefix_voi);
  voifile=fullfile(pwd,VOI_dir,'tmp_combined.voi');
else
  voifile=voifiles{1};
end
fprintf('\n');

% load VOI
voi=BVQXfile(voifile);
org_VMRdims=voi.OriginalVMRFramingCubeDim;

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
    fprintf('!!WARNING!!: some ROI(s) you specified in VOI_to_use not found. check input variable.\n');
  end
else
  voilist=1:voi.NrOfVOIs;
end

% display messages
fprintf('target VOI dir : %s\n',fullfile(pwd,VOI_dir));

[voipath,voiname,voiext]=fileparts(voifile);
fprintf('target VOI file: %s%s\n',voiname,voiext);
fprintf('\n');

% extract VOI voxel coordinates
fprintf('extracting VOI voxels...');
voi_idx=cell(numel(voilist),1);
vv=0;
for vvvv=voilist
  vv=vv+1;
  voi_voxels=voi2col(voi.VOI(vvvv).NrOfVoxels,voi.VOI(vvvv).Voxels,...
                     msk_resolution,'off',msk_box(1,:),org_VMRdims,voi.ReferenceSpace);
  tal_dims=(msk_box(2,:)-msk_box(1,:))./msk_resolution;
  voi_voxels(voi_voxels(:,1)>tal_dims(1) | voi_voxels(:,2)>tal_dims(2) | voi_voxels(:,3)>tal_dims(3),:)=[];
  voi_idx{vv}=sub2ind(tal_dims,voi_voxels(:,1),voi_voxels(:,2),voi_voxels(:,3));
end
fprintf('done.\n');
fprintf('\n');

% processing
vv=0;
msk=cell(numel(voilist),1);
fprintf('generating mask files...\n');
for vvvv=voilist
  vv=vv+1;
  fprintf('processing: %s...',voi.VOI(vvvv).Name);

  msk{vv}=BVQXfile('new:msk');
  msk{vv}.Resolution=msk_resolution;
  msk{vv}.XStart=msk_box(1,1);
  msk{vv}.XEnd=msk_box(2,1);
  msk{vv}.YStart=msk_box(1,2);
  msk{vv}.YEnd=msk_box(2,2);
  msk{vv}.ZStart=msk_box(1,3);
  msk{vv}.ZEnd=msk_box(2,3);
  msk{vv}.Mask=zeros((msk_box(2,:)-msk_box(1,:))./msk_resolution);
  msk{vv}.Mask(voi_idx{vv})=1;

  if save_flg
    fprintf('  saving...');
    msk{vv}.SaveAs(fullfile(voipath,[voi.VOI(vvvv).Name,'.msk']));
  end
  fprintf('done.\n');
end

% finally, generating an all-VOI-combined msk file
fprintf('generating all-VOI-combined mask...');
msk{vv+1}=BVQXfile('new:msk');
msk{vv+1}.Resolution=msk_resolution;
msk{vv+1}.XStart=msk_box(1,1);
msk{vv+1}.XEnd=msk_box(2,1);
msk{vv+1}.YStart=msk_box(1,2);
msk{vv+1}.YEnd=msk_box(2,2);
msk{vv+1}.ZStart=msk_box(1,3);
msk{vv+1}.ZEnd=msk_box(2,3);
msk{vv+1}.Mask=zeros((msk_box(2,:)-msk_box(1,:))./msk_resolution);
for ii=1:1:length(voi_idx), msk{vv+1}.Mask(voi_idx{ii})=1; end
msk{vv+1}.SaveAs(fullfile(voipath,'all_VOIs_combined.msk'));
fprintf('done.\n');

% clean up
voi.ClearObject(); clear voi;

if ~nargout
  for vv=1:1:length(msk), msk{vv}.ClearObject(); end
  clear msk;
end

% clean up the temporal file
if exist(fullfile(pwd,VOI_dir,'tmp_combined.voi'),'file')
  delete(fullfile(pwd,VOI_dir,'tmp_combined.voi'));
end

return

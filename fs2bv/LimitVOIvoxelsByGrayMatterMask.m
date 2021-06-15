function LimitVOIvoxelsByGrayMatterMask(VOI_dir,MSK_dir,unique_flg,prefix_voi,prefix_msk,overwrite_flg)

% Limits voxels in each VOI based on a subject-specific individual gray matter mask (required when we apply a template or group VOIs for individual analsis).
% LimitVOIvoxelsByGrayMatterMask(VOI_dir,MSK_dir,:unique_flg,:prefix_voi,:prefix_msk,:overwrite_flg)
% (: is optional)
%
% This function limits the voxels in VOI file(s) as those within the gray matter
% regions defined in the gray matter mask. This procedure would not be required
% since VOIs are generally defined in individual anatomical volume spaces. However,
% in some cases, we have to use group-averaged probablistic or template volume
% VOIs for individual analyses. In those situations, if we simply apply some
% probablistic VOIs in extracting voxel time series of a specific participant,
% the extracted voxel time series may contain responses in non-gray matter regions
% (e.g. white-matter or CSF) as anatomical incorerspondence between the template
% VOIs and a specific participant. To avoid this issue, the template VOIs should
% be confined to the individual gray matter regions in advance using this function.
%
% [input]
% VOI_dir    : Target directory that contains VOI files
%              e.g. '\CD\zk08_382\voi_files'
%              Target directory should be specified as such
%              the current directory where this function is
%              called is the origin.
% MSK_dir    : Target directory that contains MSK file to be applied.
%              e.g. '\CD\zk08_382\ROI_vtc\msk'
%              Target directory should be specified as such
%              the current directory where this function is
%              called is the origin.
% unique_flg : whether using only unique voxel values or not, 1 by default
%              if 0, the same voxel value may be counted more than once
%              1 by default.
% prefix_voi : (optional) string to determine the target VOI
%              from multiple files, e.g. '*combined_all'. empty by default.
% prefix_msk : (optional) string to determine the target MSK
%              from multiple files, e.g. 'zk08_382'. empty by default.
% overwrite_flg: (optional) whether overwrite the original file with the new one.
%                if set to 1, all the files will be overwritten by new ones.
%                if set to 0, the new data will be saved with '_limit_graymatter'
%                prefix in the same directory with the input VOIs. 0 by default.
%
% [output]
% no output variable,
% the thresholded VOIs overwrite the existing files (when overwrite_flg is 1)
% or are saved as (voi_name)_limit_graymatter.voi (when overwrite_flg is 0)
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
%
% Created    : "2012-03-15 12:10:33 banh"
% Last Update: "2018-09-03 11:32:52 ban"

% check the input variables
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(unique_flg), unique_flg=1; end
if nargin<4 || isempty(prefix_voi), prefix_voi=''; end
if nargin<5 || isempty(prefix_msk), prefix_msk=''; end
if nargin<6 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,VOI_dir),'dir')
  error('VOI_dir not found. check the input variable.');
end

if ~exist(fullfile(pwd,MSK_dir),'dir')
  error('MSK_dir not found. check the input variable.');
end

% get VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if isempty(voifiles)
  error('no voi file found. check input variable');
end

% get MSK files
mskfile=GetFiles(fullfile(pwd,MSK_dir),'*.msk',prefix_msk);
if isempty(mskfile)
  error('no msk file found. check input variable');
end
if length(mskfile)>1
  fprintf('warning: multiple MSK files are found. using the first one...\n');
end

% load MSK file and store the voxel information
[mskpath,mskname,mskext]=fileparts(mskfile{1}); %#ok
msk=BVQXfile(mskfile{1});
vox_resolution=msk.Resolution;
Xs=msk.XStart;
Ys=msk.YStart;
Zs=msk.ZStart;
%Xe=msk.XEnd;
%Ye=msk.YEnd;
%Ze=msk.ZEnd;
mask=msk.Mask;
vox_dims=size(mask);

% display message
fprintf('target VOI dir : %s\n',fullfile(pwd,VOI_dir));
fprintf('target MSK file: %s%s\n',mskname,mskext);
fprintf('\n');

% processing
for ii=1:1:length(voifiles)

  % display message
  [voipath,voifile,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s\n',voifile,voiext);

  % load VOI
  voi=BVQXfile(voifiles{ii});

  if voi.OriginalVMRResolutionX~=vox_resolution
    error('VOI and GrayMatterMask voxel resolution along X-axis mismatched. check the data structure.');
  end
  if voi.OriginalVMRResolutionY~=vox_resolution
    error('VOI and GrayMatterMask voxel resolution along Y-axis mismatched. check the data structure.');
  end
  if voi.OriginalVMRResolutionZ~=vox_resolution
    error('VOI and GrayMatterMask voxel resolution along Z-axis mismatched. check the data structure.');
  end

  org_VMRdims=voi.OriginalVMRFramingCubeDim;

  % extract VOI voxel coordinates
  fprintf('extracting VOI voxels...');
  voi_idx=cell(voi.NrOfVOIs,1);
  for vv=1:1:voi.NrOfVOIs
    voi_voxels=voi2col(voi.VOI(vv).NrOfVoxels,voi.VOI(vv).Voxels,vox_resolution,'off',[Xs,Ys,Zs],vox_dims,org_VMRdims,voi.ReferenceSpace);
    voi_voxels(voi_voxels<=0)=1; % edge adjustment, not ideal but required.
    voi_idx{vv}=sub2ind(vox_dims,voi_voxels(:,1),voi_voxels(:,2),voi_voxels(:,3));
  end
  fprintf('done.\n');

  % limit VOI voxels by the gray-matter mask
  fprintf('masking VOI voxels...');
  voxels=cell(voi.NrOfVOIs,1);
  for vv=1:1:voi.NrOfVOIs

    cvox_idx=voi_idx{vv}(logical(mask(voi_idx{vv})));
    if unique_flg, cvox_idx=unique(cvox_idx); end

    % store (x,y,z) coordinates of voxels in the whole Tal space
    [voxels{vv}(:,1),voxels{vv}(:,2),voxels{vv}(:,3)]=ind2sub(vox_dims,cvox_idx);
  end
  fprintf('done.\n');

  % updating thresholded voxels to VOI structures
  fprintf('updating voxel information...');

  if vox_resolution==1 || vox_resolution==2
    offset=[1,1,1];
  elseif vox_resolution==3
    offset=[2,2,2];
  else
    error('MSK''s voxel resolution should be one of 1,2,3.');
  end

  for vv=1:1:voi.NrOfVOIs
    voi.VOI(vv).NrOfVoxels=size(voxels{vv},1);
    tmp_vox=ones(size(voxels{vv},1),1)*ceil( 128-([Xs,Ys,Zs]-offset) )-voxels{vv}.*vox_resolution;
    voi.VOI(vv).Voxels=tmp_vox(:,[3,1,2]);
  end
  fprintf('done.\n');

  fprintf('saving the new voi...');
  if overwrite_flg
    voi.SaveAs(fullfile(voipath,[voifile,voiext]));
  else
    voi.SaveAs(fullfile(voipath,[voifile,'_limit_graymatter',voiext]));
  end
  fprintf('done.\n');

  % clean up
  voi.ClearObject(); clear voi;

end % for ii=1:1:length(voifiles)

% clean up
msk.ClearObject(); clear msk;

return

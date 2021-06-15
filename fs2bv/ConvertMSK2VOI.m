function voi=ConvertMSK2VOI(MSK_dir,voi_box,prefix_msk,save_flg)

% Converts BrainVoyaer *.msk file(s) to *.voi file(s).
% function ConvertMSK2VOI(MSK_dir,:voi_box,:prefix_msk,:save_flg)
% (: is optional)
%
% This function converts cortex-mask file(s) to a whole or partical brain VOI files.
% The generated VOI file(s) can be used for whole brain GLM, and searchlight analysis, etc.
%
% [example]
% >> voi=ConvertMSK2VOI('../fMRI_preprocessed/HB/hb14_061/ROI_vtc/mask',[256,256,256],'hb14*',1);
%
% [input]
% MSK_dir    : Target directory that contains MSK file(s)
%              e.g. '/HB/hb14_061/ROI_vtc/mask'
%              Target directory should be specified as a relative path format in which
%              the current directory where this function is called is the origin.
% voi_box    : (optional) a [1x3(x,y,z)] matrix, which define VOI box size (voxels) in the TAL space
%              if empty, [256,256,256] will be used as it is the original default MSK/VOI size in BrainVoyaer.
% prefix_msk : (optional) prefix of MSK file(s) to be processed, e.g. prefix_msk='hb14*'; empty by default.
% save_flg   : (optional) whether saving the generated voi file(s). if 1, the voi file(s)
%              wil be saved like ~/(mask_file_directory)/(mask_file_name).voi
%              0 by default.
%
% [output]
% voi        : a cell structure, voi{number_of_msk_file}. BrainVoyaer VOI object
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
% Created    : "2014-09-22 17:35:53 ban"
% Last Update: "2017-09-04 18:38:07 ban"

% check input variables
if nargin<1, help(mfilename()); if nargout, voi=[]; end; return; end
if nargin<2 || isempty(voi_box), voi_box=[256,256,256]; end
if nargin<3 || isempty(prefix_msk), prefix_msk=''; end
if nargin<4 || isempty(save_flg), save_flg=0; end

if numel(voi_box)~=3
  error('voi_box should be [1x3(x,y,z)] matrix. check input variables.');
end

if ~exist(fullfile(pwd,MSK_dir),'dir')
  error('can not find MSK_dir. check input variables.');
end

% get mask files
mskfiles=GetFiles(fullfile(pwd,MSK_dir),'*.msk',prefix_msk);
if isempty(mskfiles)
  error('can not find any MSK file in MSK_dir. check input variables.');
end

% display message
fprintf('Converting MSK file(s) to VOI(s)...\n');
fprintf('Target : %s\n\n',fullfile(pwd,MSK_dir));

% processing
voi=cell(length(mskfiles),1);
for ii=1:1:length(mskfiles)
  [mskpath,mskfname,mskext]=fileparts(mskfiles{ii});
  fprintf('processing : %s%s\n',mskfname,mskext);

  msk=BVQXfile(mskfiles{ii});
  resolution=msk.Resolution;
  if resolution==1
    offset=[1,1,1];
  elseif resolution==2
    offset=[2,2,2];
  elseif resolution==3
    offset=[2,2,2];
  else
    error('supported MSK resolution is one of 1, 2, or 3. check the input MSK: %s%s.',mskfname,mskext);
  end

  % convert the original MSK space to 1x1x1 voi space by duplicating the mask volumes along x-, y-, and z-axes.
  if resolution==1
    msk_voxels=msk.Mask;
  else
    sz=size(msk.Mask); % get the original mask box size
    msk_voxels=reshape(repmat(msk.Mask,[resolution,1,1]),[sz(1),resolution*sz(2),sz(3)]);
    msk_voxels=permute(reshape(repmat(permute(msk_voxels,[2,1,3]),[resolution,1,1]),[resolution*sz(2),resolution*sz(1),sz(3)]),[2,1,3]);
    msk_voxels=reshape(repmat(msk_voxels,[1,resolution,1]),sz*resolution);
  end

  % get voxel XYZ coordinates
  voxels=zeros(voi_box);
  voxels((msk.XStart+1:min(msk.XEnd,msk.XStart+size(msk_voxels,1)))-offset(1),...
         (msk.YStart+1:min(msk.YEnd,msk.YStart+size(msk_voxels,2)))-offset(2),...
         (msk.ZStart+1:min(msk.ZEnd,msk.ZStart+size(msk_voxels,3)))-offset(3))=msk_voxels;
  vox_idx=find(voxels~=0);
  tal_voxels=zeros(numel(vox_idx),3);
  [tal_voxels(:,2),tal_voxels(:,3),tal_voxels(:,1)]=ind2sub(voi_box,vox_idx);
  % need to caclculate 128-coordinates since the default VOI coordinates are from -128 to 128 and filpped along x-, y-, and z- axes.
  tal_voxels=128-tal_voxels;

  % generate a VOI
  voi{ii}=BVQXfile('new:voi');
  voi{ii}.VOI(1).Name=[mskfname,mskext];%strrep([mskfname,mskext],'.','_');
  voi{ii}.VOI(1).Color=[0,255,255];
  voi{ii}.VOI(1).NrOfVoxels=sum(msk_voxels);
  voi{ii}.VOI(1).Voxels=tal_voxels;

  % saving the generated voi
  if save_flg
    fprintf('saving the generate voi as %s.voi...',mskfname);
    voi{ii}.SaveAs(fullfile(mskpath,[mskfname,'.voi']));
    fprintf('done.\n');
  end

  % clean up
  msk.ClearObject(); clear msk;
  if ~nargout, voi{ii}.ClearObject(); end
  fprintf('\n');
end

return

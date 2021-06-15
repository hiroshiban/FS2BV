function FlipLeftRightMSKs(MSK_dir,prefix_msk,overwrite_flg,framingcube)

% Flips the left and right of the input BrainVoyager MSKs.
% function FlipLeftRightMSKs(MSK_dir,:prefix_msk,:overwrite_flg,:framingcube)
% (: is optional)
%
% This function loads the target MSK (mask) files and flips them along the left anf right axis.
%
% [important note]
% This function assumes that the original VMR is defined in a cubic voliume space.
% If you are using the other anisotropic spaces, please be careful.
%
% [input]
% MSK_dir    : target directory that contains MSK files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% prefix_msk : (optional) string to determine the target MSK
%              from multiple files, e.g. '*_final_TAL'. empty by default
% overwrite_flg : (optional) whether overwriting the input MSKs.
%              if 1, the original MSKs are overwritten by the flipped MSKs.
%              if 0, the flipped MSKs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flipLR'. 0 by default.
% framingcube : (optional) the volume size of the original VMRs (should be cubic).
%              256 by default. You can set 384 or 512 for high-resolution
%              scans.
%
% [output]
% no output variable.
% the flipped MSKs are saved in the original directory of the input MSKs
%
% [note on how to set the 'prefix_*' variable]
% prefix_* can be set flexibly as below.
% 1. a string: setting an including prefix (string) alone
%    e.g. prefix_*='_TDTS6.0';4
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
% Created    : "2017-12-20 14:31:10 ban"
% Last Update: "2018-09-03 12:33:35 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_msk), prefix_msk=''; end
if nargin<3 || isempty(overwrite_flg), overwrite_flg=0; end
if nargin<4 || isempty(framingcube), framingcube=256; end

if ~exist(fullfile(pwd,MSK_dir),'dir'), error('MSK_dir not found. check the input variable.'); end

% get MSK files
mskfiles=GetFiles(fullfile(pwd,MSK_dir),'*.msk',prefix_msk);
if isempty(mskfiles)
  fprintf('WARNING: no MSK file found. check the input variable. exiting the program...\n');
  return
end

% processing
for ii=1:1:length(mskfiles)
  [mskpath,mskfname,mskext]=fileparts(mskfiles{ii});
  fprintf('processing: %s%s...',mskfname,mskext);
  savefname=mskfname;
  %if ~overwrite_flg, savefname=[mskfname,'_flipZ']; end % note: flip along z-axis = flip LR
  if ~overwrite_flg, savefname=[mskfname,'_flipLR']; end
  msk=BVQXfile(mskfiles{ii});

  % flipping the coordinates
  try
    msk.Mask=flipdim(msk.Mask,3); %#ok % flip along z-axis = flip LR
  catch
    msk.Mask=flip(msk.Mask,3);
  end

  % adjust the start/end positions
  ze=framingcube-msk.ZStart;
  zs=framingcube-msk.ZEnd;
  msk.ZEnd=ze;
  msk.ZStart=zs;

  fprintf('done.\n');

  fprintf('saving    : %s%s...',savefname,mskext);
  msk.SaveAs(fullfile(mskpath,[savefname,mskext]));
  fprintf('done.\n');

  msk.ClearObject(); clear msk;
end

return

function FlipLeftRightSRFs(SRF_dir,prefix_srf,overwrite_flg)

% Flips the left and right of the input BrainVoyager SRFs.
% function FlipLeftRightSRFs(SRF_dir,:prefix_srf,:overwrite_flg)
% (: is optional)
%
% This function is a simple wrapper of FlipSRFs to be used as
% one of a series of functions to flip LR of the BrainVoyager
% datasets (for details, please see FlipSRFs).
% Briefly, this function loads the target BrainVoyager SRF surface
% files and flips their left and right coordinates.
%
% [note]
% This function also flips the left and right of the linked SRF file
% names (_LH <--> _RH) as well as the coordinates of the surface vertices.
% Therefore, if you only flip *_LH.srf alone, some inconsistency of file
% structures may happen (it is actually not a big issue though. the linked
% reference surface file is just simply detached). Please be careful.
%
% [input]
% SRF_dir    : target directory that contains SRF files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% prefix_srf : (optional) string to determine the target SRF
%              from multiple files, e.g. '*_final_TAL'. empty by default
% overwrite_flg : (optional) whether overwriting the input SRFs.
%              if 1, the original SRFs are overwritten by the flipped SRFs.
%              if 0, the flipped SRFs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flipLR'. 0 by default.
%
% [output]
% no output variable.
% the flipped SRFs are saved in the original directory of the input SRFs
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
% Created    : "2017-12-21 10:53:18 ban"
% Last Update: "2017-12-23 16:06:48 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_srf), prefix_srf=''; end
if nargin<3 || isempty(overwrite_flg), overwrite_flg=0; end

FlipSRFs(SRF_dir,[0,0,1],prefix_srf,overwrite_flg); % [0,0,1] = flip along the z-axis = flip LR;

return

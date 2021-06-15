function FlipLeftRightVMRs(VMR_dir,prefix_vmr,overwrite_flg)

% Flips the left and right of the input BrainVoyager VMRs and the correponding V16.
% function FlipLeftRightVMRs(VMR_dir,:prefix_vmr,:overwrite_flg)
% (: is optional)
%
% This function is a simple wrapper of FlipVMRs to be used as
% one of a series of functions to flip LR of the BrainVoyager
% datasets (for details, please see FlipVMRs).
% Briefly, this function loads the target BrainVoyager VMR volume
% files and flips their left and right coordinates.
%
% [input]
% VMR_dir    : target directory that contains VMR files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% prefix_vmr : (optional) string to determine the target VMR
%              from multiple files, e.g. '*_final_TAL'. empty by default
% overwrite_flg : (optional) whether overwriting the input VMRs.
%              if 1, the original VMRs are overwritten by the flipped VMRs.
%              if 0, the flipped VMRs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flipLR'. 0 by default.
%
% [output]
% no output variable.
% the flipped VMRs are saved in the original directory of the input VMRs
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
% Created    : "2014-03-24 14:00:45 ban"
% Last Update: "2017-12-23 15:19:47 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_vmr), prefix_vmr=''; end
if nargin<3 || isempty(overwrite_flg), overwrite_flg=0; end

FlipVMRs(VMR_dir,[0,0,1],prefix_vmr,overwrite_flg); % [0,0,1] = flip along the z-axis = flip LR

return

function ShiftVOIs(VOI_dir,shift_XYZ,prefix_voi,overwrite_flg)

% Translates (parallel transport) VOI coordinates along XYZ-axes.
% function ShiftVOIs(VOI_dir,shift_XYZ,:prefix_voi,:overwrite_flg)
% (: is optional)
%
% This function shifts VOI positions parallelly along XYZ-axes.
%
% [input]
% VOI_dir    : Target directory that contains VOI files
%              e.g. '\HB\hb17_010\voi_files'
%              Target directory should be specified as a relative path
%              format in which the current directory where this function is
%              called is the origin.
% shift_XYZ  : amounts of shift along XYZ, [X,Y,Z], which should be set
%              by voxel base (1 = 1 voxel shift). e.g. shift_XYZ=[-3,0,0];
% prefix_voi : (optional) string to determine the target from
%              multiple files, e.g. '_combined'.
%              empty by default.
% overwrite_flg: (optional) whether overwrite the original file with the new one.
%              if set to 1, all the files will be overwritten by new ones.
%              if set to 0, the new data will be saved with '_vmpmask' prefix
%              in the same directory. 0 by default.
%
% [output]
% no output variable
% The translated VOI will overwrite the existing one (when overwrite_flg is 1)
% or is saved with "_shifted" prefix (when overwrite_flg is 0)
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
% 2. wildcardsearch.m
% enable reg-exp search of files
%
%
% Created    : "2021-05-28 17:38:42 ban"
% Last Update: "2021-05-28 18:15:07 ban"

% check the input variable
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(prefix_voi), prefix_voi=''; end
if nargin<4 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,VOI_dir),'dir')
  error('VOI_dir not found. check the input variable.');
end

if numel(shift_XYZ)~=3
  error('shift_XYZ should be [X,Y,Z]. check the input variable.');
end

fprintf('Translating VOI by [% 3d, % 3d, %3d]...\n',shift_XYZ(1),shift_XYZ(2),shift_XYZ(3));
fprintf('Traget : %s\n',VOI_dir);

% get the target VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);

% processing
for ii=1:1:length(voifiles)
  [voipath,voifname,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s...',voifname,voiext);

  % load the target VOI file
  voi=BVQXfile(voifiles{ii});

  % translating
  for rr=1:1:voi.NrOfVOIs
    voi.VOI(rr).Voxels=voi.VOI(rr).Voxels+repmat([shift_XYZ(1),shift_XYZ(2),shift_XYZ(3)],[size(voi.VOI(rr).Voxels,1),1]);
  end

  % saving
  tgtvoipath=relativepath(voipath);
  tgtvoipath=tgtvoipath(1:end-1); %#ok
  if overwrite_flg
    voi.SaveAs(fullfile(voipath,[voifname,voiext]));
    evalc('CompactVOIs(tgtvoipath,{{voifname};{''_shifted''}})');
  else
    voi.SaveAs(fullfile(voipath,[voifname,'_shifted',voiext]));
    evalc('CompactVOIs(tgtvoipath,[voifname,''_shifted''])');
  end

  % clean up
  voi.ClearObject();
  clear voi;

  fprintf('done.\n');
end

fprintf('completed.\n');

return

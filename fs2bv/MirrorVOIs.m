function MirrorVOIs(VOI_dir,prefix_voi)

% Mirros VOI. Specifically, the VOI corrdinates are flipped along x-axis.
% function MirrorVOIs(VOI_dir,:prefix_voi)
% (: is optional)
%
% This function mirrors VOIs by flipping the voxel coordinates along left/right-axis.
%
% [note on the difference between MirrorVOIs and FlipLeftRightVOIs]
% MirrorVOIs        : just flip the left and right without changing the VOI names etc
% FlipLeftRightVOIs : this flips all the data structures, including to change '_lh' and
%                     '_rh' of the internally stored VOI names.
%
% [input]
% VOI_dir    : Target directory that contains VOI files
%              e.g. '\HB\hb17_010\voi_files'
%              Target directory should be specified as a relative path
%              format in which the current directory where this function is
%              called is the origin.
% prefix_voi : (optional) string to determine the target from
%              multiple files, e.g. '_combined'.
%              empty by default.
%
% [output]
% no output variable
% the flipped VOI files are saved with a '_mirrored' prefix, as (voi_file_name)_mirrored.voi.
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
% Created    : "2017-09-07 09:55:24 ban"
% Last Update: "2017-12-21 10:57:43 ban"

% check the input variable
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_voi), prefix_voi=''; end

if ~exist(fullfile(pwd,VOI_dir),'dir')
  error('VOI_dir not found. check the input variable.');
end

fprintf('Mirroring voxels (flipping coordinates along x-axis)...\n');
fprintf('Traget : %s\n',VOI_dir);

% get the target VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);

% processing
for ii=1:1:length(voifiles)
  [voipath,voifname,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s...',voifname,voiext);

  % load the target VOI file
  voi=BVQXfile(voifiles{ii});

  % mirroring the voxels
  for vv=1:1:voi.NrOfVOIs
    voi.VOI(vv).Voxels(:,1)=-1*voi.VOI(vv).Voxels(:,1);
  end

  % saving
  %voi.SaveAs(fullfile(voipath,[voifname,'_mirrored',voiext]));

  % updated
  % fopen the same target VOI file to overwrite
  fid=fopen(voifiles{ii},'W');
  if fid==-1, error('%s%s can not be open. check the file, run fclose all, and try again.',viofname,voiext); end

  % write the VOI details into the file
  fprintf(fid,'\nFileVersion:                %d\n',voi.FileVersion);

  if voi.FileVersion<4
    fprintf(fid,'\nCoordsType:                 %s\n',voi.ReferenceSpace);
    fprintf(fid,'\nSubjectVOINamingConvention: %s\n',voi.SubjectVOINamingConvention);
    fprintf(fid,'\n');
  else
    fprintf(fid,'\nRefereneSpace:              %s\n',voi.ReferenceSpace);
    fprintf(fid,'\n');
    fprintf(fid,'OriginalVMRResolutionX:     %d\n',voi.OriginalVMRResolutionX);
    fprintf(fid,'OriginalVMRResolutionY:     %d\n',voi.OriginalVMRResolutionY);
    fprintf(fid,'OriginalVMRResolutionZ:     %d\n',voi.OriginalVMRResolutionZ);
    fprintf(fid,'OriginalVMROffsetX:         %d\n',voi.OriginalVMROffsetX);
    fprintf(fid,'OriginalVMROffsetY:         %d\n',voi.OriginalVMROffsetY);
    fprintf(fid,'OriginalVMROffsetZ:         %d\n',voi.OriginalVMROffsetZ);
    fprintf(fid,'OriginalVMRFramingCubeDim:  %d\n',voi.OriginalVMRFramingCubeDim);
    fprintf(fid,'\nLeftRightConvention:        %d\n',voi.Convention);
    fprintf(fid,'\nSubjectVOINamingConvention: %s\n',voi.SubjectVOINamingConvention);
  end

  fprintf(fid,'\nNrOfVOIs:                   %d\n',voi.NrOfVOIs);
  fprintf(fid,'\n');

  for vv=1:1:voi.NrOfVOIs
    fprintf(fid,'NameOfVOI:                  %s\n',voi.VOI(vv).Name);
    fprintf(fid,'ColorOfVOI:                 %d %d %d\n',voi.VOI(vv).Color(1),voi.VOI(vv).Color(2),voi.VOI(vv).Color(3));
    fprintf(fid,'\n');
    fprintf(fid,'NrOfVoxels:                 %d\n',voi.VOI(vv).NrOfVoxels);
    voxelfmt='%d %d %d\n'; % * HERE, REMOVING THE REDUNDANT WHITE SPACES *
    fprintf(fid,voxelfmt,(voi.VOI(vv).Voxels)');
    fprintf(fid,'\n\n');
  end
  fprintf(fid,'NrOfVOIVTCs:                %d',voi.NrOfVTCs);

  % fclose the VOI file
  fclose(fid);

  % clean up
  voi.ClearObject(); clear voi;

  fprintf('done.\n');
end

fprintf('completed.\n');

return

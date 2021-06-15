function FlipLeftRightVOIs(VOI_dir,prefix_voi,overwrite_flg)

% Flips the left and right of the input BrainVoyager VOIs.
% function FlipLeftRightVOIs(VOI_dir,:prefix_voi,:overwrite_flg)
% (: is optional)
%
% This function loads the target VOI (Volumes Of Interest = ROI) files and flips them along the left anf right axis.
%
% [note]
% This function also flips the left and right of the VOI names as well as the VOI left/right coordinates.
% For instance, a VOI defined as V1_lh is renamed to V1_rh. Please be careful.
%
% [note on the difference between MirrorVOIs and FlipLeftRightVOIs]
% MirrorVOIs        : just flip the left and right without changing the VOI names etc
% FlipLeftRightVOIs : this flips all the data structure, including to change '_lh' and
%                     '_rh' of the internally stored VOI names.
%
% [input]
% VOI_dir    : target directory that contains VOI files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% prefix_voi : (optional) string to determine the target VOI
%              from multiple files, e.g. 'HB_hb17_010*'. empty by default
% overwrite_flg : (optional) whether overwriting the input VOIs.
%              if 1, the original VOIs are overwritten by the flipped VOIs.
%              if 0, the flipped VOIs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flipLR'. 0 by default.
%
% [output]
% no output variable.
% the flipped VOIs are saved in the original directory of the input VOIs
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
% Created    : "2017-12-20 16:35:25 ban"
% Last Update: "2018-09-03 12:34:51 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_voi), prefix_voi=''; end
if nargin<3 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,VOI_dir),'dir'), error('VOI_dir not found. check the input variable.'); end

% get VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if isempty(voifiles)
  fprintf('WARNING: no VOI file found. check the input variable. exiting the program...\n');
  return
end

% processing
for ii=1:1:length(voifiles)
  [voipath,voifname,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s...',voifname,voiext);
  savefname=voifname;
  %if ~overwrite_flg, savefname=[voifname,'_flipZ']; end % note: flip along z-axis = flip LR
  if ~overwrite_flg, savefname=[voifname,'_flipLR']; end
  voi=BVQXfile(voifiles{ii});

  for vv=1:1:voi.NrOfVOIs
    % flipping the coordinates
    if strcmpi(voi.ReferenceSpace,'TAL')
      voi.VOI(vv).Voxels(:,1)=-1*voi.VOI(vv).Voxels(:,1);
    else
      voi.VOI(vv).Voxels(:,1)=voi.OriginalVMRFramingCubeDim-voi.VOI(vv).Voxels(:,1);
    end

    % flipping the VOI names
    if strfind(voi.VOI(vv).Name,'_lh')
      voi.VOI(vv).Name=strrep(voi.VOI(vv).Name,'_lh','_rh');
    elseif strfind(voi.VOI(vv).Name,'_LH')
      voi.VOI(vv).Name=strrep(voi.VOI(vv).Name,'_LH','_RH');
    elseif strfind(voi.VOI(vv).Name,'_rh')
      voi.VOI(vv).Name=strrep(voi.VOI(vv).Name,'_rh','_lh');
    elseif strfind(voi.VOI(vv).Name,'_RH')
      voi.VOI(vv).Name=strrep(voi.VOI(vv).Name,'_RH','_LH');
    end
  end

  fprintf('done.\n');

  % display warning
  %fprintf('*** the VOI name was also flipped [lh <--> rh] as well as the XYZ coordinates. please be careful. ***\n');

  fprintf('saving    : %s%s...',savefname,voiext);
  if ~overwrite_flg
    fid=fopen(fullfile(voipath,[savefname,voiext]),'W');
  else
    % swap the file prefix '_LH' and '_RH'.
    % Here, the files are renamed with new prefixes '_lhlh' and '_rhrh' to
    % avoid overwriting the existing files to be processed in the loop.
    if ~isempty(strfind(savefname,'_lh'))
      fid=fopen(fullfile(voipath,[strrep(savefname,'_lh','_rhrh'),voiext]),'W');
    elseif ~isempty(strfind(savefname,'_LH'))
      fid=fopen(fullfile(voipath,[strrep(savefname,'_LH','_RHRH'),voiext]),'W');
    elseif ~isempty(strfind(savefname,'_rh'))
      fid=fopen(fullfile(voipath,[strrep(savefname,'_rh','_lhlh'),voiext]),'W');
    elseif ~isempty(strfind(savefname,'_RH'))
      fid=fopen(fullfile(voipath,[strrep(savefname,'_RH','_LHLH'),voiext]),'W');
    else % just overwite
      fid=fopen(fullfile(voipath,[savefname,voiext]),'W');
    end
  end
  if fid==-1, error('%s%s can not be open. check the file, run fclose all, and try again.',voifname,voiext); end

  % write the VOI details into the file
  fprintf(fid,'\nFileVersion:                %d\n',voi.FileVersion);

  if voi.FileVersion<4
    fprintf(fid,'\nCoordsType:                 %s\n',voi.ReferenceSpace);
    fprintf(fid,'\nSubjectVOINamingConvention: %s\n',voi.SubjectVOINamingConvention);
    fprintf(fid,'\n');
  else
    fprintf(fid,'\nReferenceSpace:             %s\n',voi.ReferenceSpace);
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

  fprintf('done.\n');

  voi.ClearObject(); clear voi;
end

% renaming the processed files with the proper left and right prefixes
if overwrite_flg
  for ii=1:1:length(voifiles)
    [voipath,voifname,voiext]=fileparts(voifiles{ii});
    if ~isempty(strfind(voifname,'_lh'))
      movefile(fullfile(voipath,[strrep(voifname,'_lh','_rhrh'),voiext]),fullfile(voipath,[strrep(strrep(voifname,'_lh','_rhrh'),'_rhrh','_rh'),voiext]),'f');
    elseif ~isempty(strfind(voifname,'_LH'))
      movefile(fullfile(voipath,[strrep(voifname,'_LH','_RHRH'),voiext]),fullfile(voipath,[strrep(strrep(voifname,'_LH','_RHRH'),'_RHRH','_RH'),voiext]),'f');
    elseif ~isempty(strfind(voifname,'_rh'))
      movefile(fullfile(voipath,[strrep(voifname,'_rh','_lhlh'),voiext]),fullfile(voipath,[strrep(strrep(voifname,'_rh','_lhlh'),'_lhlh','_lh'),voiext]),'f');
    elseif ~isempty(strfind(voifname,'_RH'))
      movefile(fullfile(voipath,[strrep(voifname,'_RH','_LHLH'),voiext]),fullfile(voipath,[strrep(strrep(voifname,'_RH','_LHLH'),'_LHLH','_LH'),voiext]),'f');
    end
  end
end

return

function FlipSRFs(SRF_dir,XYZ_flg,prefix_srf,overwrite_flg)

% Flips the input BrainVoyager SRF surfaces along X/Y/Z-axis.
% function FlipLeftRightSRFs(SRF_dir,:XYZ_flg,:prefix_srf,:overwrite_flg)
% (: is optional)
%
% This function loads the target BrainVoyager SRF surface files and flips them along X/Y/Z-axis.
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
% XYZ_flg    : (optional) flags to specify the axis(axes) along which
%              the images are flipped. 1x3(XYZ) vector. 1 means to be
%              flipped along the corresponding axis.
%              XYZ_flg=[0,0,1]; by default.
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
% Created    : "2017-12-21 10:59:08 ban"
% Last Update: "2018-09-03 12:35:09 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(XYZ_flg), XYZ_flg=[0,0,1]; end
if nargin<3 || isempty(prefix_srf), prefix_srf=''; end
if nargin<4 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,SRF_dir),'dir'), error('SRF_dir not found. check the input variable.'); end

% get SRF files
srffiles=GetFiles(fullfile(pwd,SRF_dir),'*.srf',prefix_srf);
if isempty(srffiles)
  fprintf('WARNING: no SRF file found. check the input variable. exiting the program...\n');
  return
end

% processing
savefname=cell(length(srffiles),1);
for ii=1:1:length(srffiles)
  [srfpath,srffname,srfext]=fileparts(srffiles{ii});
  fprintf('processing: %s%s...',srffname,srfext);
  savefname{ii}=srffname;
  if ~overwrite_flg, savefname{ii}=[srffname,'_flip']; end
  srf=BVQXfile(srffiles{ii});

  % flipping the coordinates
  if XYZ_flg(1)
    srf=srf.Transform([-1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]); % affine transformation, flip the 3D matrix along the x-axis
    if ~overwrite_flg, savefname{ii}=[savefname{ii},'X']; end
  end

  if XYZ_flg(2)
    srf=srf.Transform([1,0,0,0;0,-1,0,0;0,0,1,0;0,0,0,1]); % affine transformation, flip the 3D matrix along the y-axis
    if ~overwrite_flg, savefname{ii}=[savefname{ii},'Y']; end
  end

  if XYZ_flg(3)
    srf=srf.Transform([1,0,0,0;0,1,0,0;0,0,-1,0;0,0,0,1]); % affine transformation, flip the 3D matrix along the z-axis (= LR)
    if ~overwrite_flg, savefname{ii}=[savefname{ii},'Z']; end
  end

  srf.TriangleVertex=srf.TriangleVertex(:,[3,2,1]);
  warning off; %#ok
  srf.Neighbors=srf.TrianglesToNeighbors(srf.VertexCoordinate);
  warning on; %#ok
  srf=srf.RecalcNormals();

  fprintf('done.\n');

  % display warning on the reference recosm surface (srf.AutoLinkedMTC) change.
  % the actual updates of the reference SRF is done later
  %srf.AutoLinkedMTC='';
  if ~isempty(srf.AutoLinkedMTC)
    fprintf('*** the reference surface name is also flipped [lh <--> rh]. please be careful. ***\n');
  end

  % setup for a special flipping case -- left/right flip
  if XYZ_flg(1)==0 && XYZ_flg(2)==0 && XYZ_flg(3)==1
    if ~overwrite_flg, savefname{ii}=[srffname,'_flipLR']; end
  end

  fprintf('saving    : %s%s...',savefname{ii},srfext);
  if ~overwrite_flg
    srf.SaveAs(fullfile(srfpath,[savefname{ii},srfext]));
  else
    % swap the file prefix '_LH' and '_RH'.
    % Here, the files are renamed with new prefixes '_lhlh' and '_rhrh' to
    % avoid overwriting the existing files to be processed in the loop.
    if ~isempty(strfind(savefname{ii},'_lh'))
      srf.SaveAs(fullfile(srfpath,[strrep(savefname{ii},'_lh','_rhrh'),srfext]));
    elseif ~isempty(strfind(savefname{ii},'_LH'))
      srf.SaveAs(fullfile(srfpath,[strrep(savefname{ii},'_LH','_RHRH'),srfext]));
    elseif ~isempty(strfind(savefname{ii},'_rh'))
      srf.SaveAs(fullfile(srfpath,[strrep(savefname{ii},'_rh','_lhlh'),srfext]));
    elseif ~isempty(strfind(savefname{ii},'_RH'))
      srf.SaveAs(fullfile(srfpath,[strrep(savefname{ii},'_RH','_LHLH'),srfext]));
    else % just overwite
      srf.SaveAs(fullfile(srfpath,[savefname{ii},srfext]));
    end
  end
  fprintf('done.\n');

  srf.ClearObject(); clear srf;
end

% renaming the processed files with the proper left and right prefixes
if overwrite_flg
  for ii=1:1:length(srffiles)
    [srfpath,srffname,srfext]=fileparts(srffiles{ii});
    if ~isempty(strfind(srffname,'_lh'))
      movefile(fullfile(srfpath,[strrep(srffname,'_lh','_rhrh'),srfext]),fullfile(srfpath,[strrep(strrep(srffname,'_lh','_rhrh'),'_rhrh','_rh'),srfext]),'f');
    elseif ~isempty(strfind(srffname,'_LH'))
      movefile(fullfile(srfpath,[strrep(srffname,'_LH','_RHRH'),srfext]),fullfile(srfpath,[strrep(strrep(srffname,'_LH','_RHRH'),'_RHRH','_RH'),srfext]),'f');
    elseif ~isempty(strfind(srffname,'_rh'))
      movefile(fullfile(srfpath,[strrep(srffname,'_rh','_lhlh'),srfext]),fullfile(srfpath,[strrep(strrep(srffname,'_rh','_lhlh'),'_lhlh','_lh'),srfext]),'f');
    elseif ~isempty(strfind(srffname,'_RH'))
      movefile(fullfile(srfpath,[strrep(srffname,'_RH','_LHLH'),srfext]),fullfile(srfpath,[strrep(strrep(srffname,'_RH','_LHLH'),'_LHLH','_LH'),srfext]),'f');
    end
  end
end

% re-attaching the reference SRF files
if overwrite_flg
  tmpsrf=GetFiles(SRF_dir,SRF_dir,'*RECOSM.srf','LH');
  if ~isempty(tmpsrf)
    AttachRECOSM2SRF(SRF_dir,SRF_dir,'LH','LH');
  end
  tmpsrf=GetFiles(SRF_dir,SRF_dir,'*RECOSM.srf','RH');
  if ~isempty(tmpsrf)
    AttachRECOSM2SRF(SRF_dir,SRF_dir,'RH','RH');
  end
else
  tmpsrf=GetFiles(SRF_dir,SRF_dir,'*_RECOSM.srf','LH*flip');
  if ~isempty(tmpsrf)
    AttachRECOSM2SRF(SRF_dir,SRF_dir,'LH*flip','LH*flip');
  end
  tmpsrf=GetFiles(SRF_dir,SRF_dir,'*_RECOSM.srf','RH*flip');
  if ~isempty(tmpsrf)
    AttachRECOSM2SRF(SRF_dir,SRF_dir,'RH*flip','RH*flip');
  end
end

return

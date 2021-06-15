function FlipLeftRightPOIs(POI_dir,prefix_poi,overwrite_flg)

% Flips the left and right of the input BrainVoyager POIs.
% function FlipLeftRightPOIs(POI_dir,:prefix_poi,:overwrite_flg)
% (: is optional)
%
% This function loads the target POI (Patches Of Interest = surface ROI) files and flips them along the left anf right axis.
%
% [note]
% Internally, flipping of the POIs are done by
% 1. changing the prefixes of the linked SRF (smp.NameOfOriginalSRF) names (_LH <--> _RH), and
% 2. changing the left and right of the POI names. For instance, a VOI defined as V1_lh is renamed to V1_rh.
% Please be careful.
%
% [input]
% POI_dir    : target directory that contains POI files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% prefix_poi : (optional) string to determine the target POI
%              from multiple files, e.g. 'HB_hb17_010*'. empty by default
% overwrite_flg : (optional) whether overwriting the input POIs.
%              if 1, the original POIs are overwritten by the flipped POIs.
%              if 0, the flipped POIs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flipLR'. 0 by default.
%
% [output]
% no output variable.
% the flipped POIs are saved in the original directory of the input POIs
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
% Created    : "2017-12-20 17:38:42 ban"
% Last Update: "2018-09-03 12:33:43 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_poi), prefix_poi=''; end
if nargin<3 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,POI_dir),'dir'), error('POI_dir not found. check the input variable.'); end

% get POI files
poifiles=GetFiles(fullfile(pwd,POI_dir),'*.poi',prefix_poi);
if isempty(poifiles)
  fprintf('WARNING: no POI file found. check the input variable. exiting the program...\n');
  return
end

% processing
for ii=1:1:length(poifiles)
  [poipath,poifname,poiext]=fileparts(poifiles{ii});
  fprintf('processing: %s%s...',poifname,poiext);
  savefname=poifname;
  %if ~overwrite_flg, savefname=[poifname,'_flipZ']; end % note: flip along z-axis = flip LR
  if ~overwrite_flg, savefname=[poifname,'_flipLR']; end
  poi=BVQXfile(poifiles{ii});

  % flipping the POI names
  for vv=1:1:poi.NrOfPOIs
    if strfind(poi.POI(vv).Name,'_lh')
      poi.POI(vv).Name=strrep(poi.POI(vv).Name,'_lh','_rh');
    elseif strfind(poi.POI(vv).Name,'_LH')
      poi.POI(vv).Name=strrep(poi.POI(vv).Name,'_LH','_RH');
    elseif strfind(poi.POI(vv).Name,'_rh')
      poi.POI(vv).Name=strrep(poi.POI(vv).Name,'_rh','_lh');
    elseif strfind(poi.POI(vv).Name,'_RH')
      poi.POI(vv).Name=strrep(poi.POI(vv).Name,'_RH','_LH');
    end
  end

  % flipping the reference
  if ~isempty(strfind(poi.FromMeshFile,'_LH'))
    poi.FromMeshFile=strrep(strrep(poi.FromMeshFile,'_LH','_RH'),filesep(),'/');
  elseif ~isempty(strfind(poi.FromMeshFile,'_lh'))
    poi.FromMeshFile=strrep(strrep(poi.FromMeshFile,'_lh','_rh'),filesep(),'/');
  elseif ~isempty(strfind(poi.FromMeshFile,'_RH'))
    poi.FromMeshFile=strrep(strrep(poi.FromMeshFile,'_RH','_LH'),filesep(),'/');
  elseif ~isempty(strfind(poi.FromMeshFile,'_rh'))
    poi.FromMeshFile=strrep(strrep(poi.FromMeshFile,'_rh','_lh'),filesep(),'/');
  end

  fprintf('done.\n');

  % display warning
  %fprintf('*** the POI name was also flipped [lh <--> rh] as well as the XYZ coordinates. please be careful. ***\n');

  fprintf('saving    : %s%s...',savefname,poiext);
  if ~overwrite_flg
    poi.SaveAs(fullfile(poipath,[savefname,poiext]));
  else
    % swap the file prefix '_LH' and '_RH'.
    % Here, the files are renamed with new prefixes '_lhlh' and '_rhrh' to
    % avoid overwriting the existing files to be processed in the loop.
    if ~isempty(strfind(poifname,'_lh'))
      poi.SaveAs(fullfile(poipath,[strrep(savefname,'_lh','_rhrh'),poiext]));
    elseif ~isempty(strfind(poifname,'_LH'))
      poi.SaveAs(fullfile(poipath,[strrep(savefname,'_LH','_RHRH'),poiext]));
    elseif ~isempty(strfind(poifname,'_rh'))
      poi.SaveAs(fullfile(poipath,[strrep(savefname,'_rh','_lhlh'),poiext]));
    elseif ~isempty(strfind(poifname,'_RH'))
      poi.SaveAs(fullfile(poipath,[strrep(savefname,'_RH','_LHLH'),poiext]));
    else % just overwite
      poi.SaveAs(fullfile(poipath,[savefname,poiext]));
    end
  end
  fprintf('done.\n');

  poi.ClearObject(); clear poi;
end

% renaming the processed files with the proper left and right prefixes
if overwrite_flg
  for ii=1:1:length(poifiles)
    [poipath,poifname,poiext]=fileparts(poifiles{ii});
    if ~isempty(strfind(poifname,'_lh'))
      movefile(fullfile(poipath,[strrep(poifname,'_lh','_rhrh'),poiext]),fullfile(poipath,[strrep(strrep(poifname,'_lh','_rhrh'),'_rhrh','_rh'),poiext]),'f');
    elseif ~isempty(strfind(poifname,'_LH'))
      movefile(fullfile(poipath,[strrep(poifname,'_LH','_RHRH'),poiext]),fullfile(poipath,[strrep(strrep(poifname,'_LH','_RHRH'),'_RHRH','_RH'),poiext]),'f');
    elseif ~isempty(strfind(poifname,'_rh'))
      movefile(fullfile(poipath,[strrep(poifname,'_rh','_lhlh'),poiext]),fullfile(poipath,[strrep(strrep(poifname,'_rh','_lhlh'),'_lhlh','_lh'),poiext]),'f');
    elseif ~isempty(strfind(poifname,'_RH'))
      movefile(fullfile(poipath,[strrep(poifname,'_RH','_LHLH'),poiext]),fullfile(poipath,[strrep(strrep(poifname,'_RH','_LHLH'),'_LHLH','_LH'),poiext]),'f');
    end
  end
end

return

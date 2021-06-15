function FlipVMRs(VMR_dir,XYZ_flg,prefix_vmr,overwrite_flg)

% Flips BrainVoyager VMRs and the corresponding V16 along X/Y/Z-axis.
% function FlipVMRs(VMR_dir,:XYZ_flg,:prefix_vmr,:overwrite_flg)
% (: is optional)
%
% This function loads the target BrainVoyager VMR volume files and
% flips them along X/Y/Z-axis.
%
% [input]
% VMR_dir    : target directory that contains VMR files
%              the directory should be specified with a relative path
%              format so that the current directory where this function
%              is called is the origin.
% XYZ_flg    : (optional) flags to specify the axis(axes) along which
%              the images are flipped. 1x3(X: Anterior-Posterior,
%              Y: Superior-Inferior, Z: Left-Right) vector. 1 means to be
%              flipped along the corresponding axis.
%              XYZ_flg=[0,0,1]; by default.
% prefix_vmr : (optional) string to determine the target VMR
%              from multiple files, e.g. '*_final_TAL'. empty by default
% overwrite_flg : (optional) whether overwriting the input VMRs.
%              if 1, the original VMRs are overwritten by the flipped VMRs.
%              if 0, the flipped VMRs will be stored in the same directory
%              with the location where the original ones are, with a file
%              prefix '_flip(X)(Y)(Z)'. 0 by default.
%
% [output]
% no output variable.
% the flipped VMRs are stored in the original directory where the input
% VMRs are stored.
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
% Last Update: "2018-09-03 12:35:14 ban"

% check input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(XYZ_flg), XYZ_flg=[0,0,1]; end
if nargin<3 || isempty(prefix_vmr), prefix_vmr=''; end
if nargin<4 || isempty(overwrite_flg), overwrite_flg=0; end

if ~exist(fullfile(pwd,VMR_dir),'dir'), error('VMR_dir not found. check the input variable.'); end

% get VMR files
vmrfiles=GetFiles(fullfile(pwd,VMR_dir),'*.vmr',prefix_vmr);
if isempty(vmrfiles)
  fprintf('WARNING: no VMR file found. check the input variable. exiting the program...\n');
  return
end

% processing
for ii=1:1:length(vmrfiles)
  [vmrpath,vmrfname,vmrext]=fileparts(vmrfiles{ii});
  fprintf('processing: %s%s...',vmrfname,vmrext);
  savefname=vmrfname;
  if ~overwrite_flg, savefname=[vmrfname,'_flip']; end
  vmr=BVQXfile(vmrfiles{ii});

  if XYZ_flg(1) % along anterior-posterior
    try
      vmr.VMRData=flipdim(vmr.VMRData,1);     %#ok
      vmr.VMRData16=flipdim(vmr.VMRData16,1); %#ok
    catch
      vmr.VMRData=flip(vmr.VMRData,1);
      vmr.VMRData16=flip(vmr.VMRData16,1);
    end

    % flip slice position
    vmr.OffsetX=-1*vmr.OffsetX;
    vmr.Slice1CenterX=-1*vmr.Slice1CenterX;
    vmr.SliceNCenterX=-1*vmr.SliceNCenterX;
    if ~overwrite_flg, savefname=[savefname,'X']; end %#ok

    % flip the transformation matrix
    for vv=1:1:numel(vmr.Trf,1)
      if vmr.Trf(vv).NrOfSpatialTransformationValues==16 && ( size(vmr.Trf(vv).TransformationValues,1)==4 && size(vmr.Trf(vv).TransformationValues,2)==4 )
        vmr.Trf(vv).TransformationValues=[1,-1,-1,-1; -1,1,1,1; -1,1,1,1; 1,1,1,1].*vmr.Trf(vv).TransformationValues;
      elseif vmr.Trf(vv).NrOfSpatialTransformationValues==24 && size(vmr.Trf(vv).TransformationValues,1)==24 % Talairach transformation
        vmr.Trf(vv).TransformationValues(1:3:24)=...
            vmr.Trf(vv).TransformationValues(1)-(vmr.Trf(vv).TransformationValues(1:3:24)-vmr.Trf(vv).TransformationValues(1));
        % swap anterior and posterior coordinates
        tmp=vmr.Trf(vv).TransformationValues(10);
        vmr.Trf(vv).TransformationValues(10)=vmr.Trf(vv).TransformationValues(10-3);
        vmr.Trf(vv).TransformationValues(10-3)=tmp;
        clear tmp;
      else
        error('the transformation matrix in %s%s is not compatible wiht this function.',vmrfname,vmrext);
      end
    end
  end

  if XYZ_flg(2) % along superior-inferior
    try
      vmr.VMRData=flipdim(vmr.VMRData,2);     %#ok
      vmr.VMRData16=flipdim(vmr.VMRData16,2); %#ok
    catch
      vmr.VMRData=flip(vmr.VMRData,2);
      vmr.VMRData16=flip(vmr.VMRData16,2);
    end

    % flip slice position
    vmr.OffsetY=-1*vmr.OffsetY;
    vmr.Slice1CenterY=-1*vmr.Slice1CenterY;
    vmr.SliceNCenterY=-1*vmr.SliceNCenterY;
    if ~overwrite_flg, savefname=[savefname,'Y']; end %#ok

    % flip the transformation matrix
    for vv=1:1:numel(vmr.Trf,1)
      if vmr.Trf(vv).NrOfSpatialTransformationValues==16 && ( size(vmr.Trf(vv).TransformationValues,1)==4 && size(vmr.Trf(vv).TransformationValues,2)==4 )
        vmr.Trf(vv).TransformationValues=[1,-1,1,1; -1,1,-1,-1; 1,-1,1,1; 1,1,1,1].*vmr.Trf(vv).TransformationValues;
      elseif vmr.Trf(vv).NrOfSpatialTransformationValues==24 && size(vmr.Trf(vv).TransformationValues,1)==24 % Talairach transformation
        vmr.Trf(vv).TransformationValues(2:3:24)=...
            vmr.Trf(vv).TransformationValues(2)-(vmr.Trf(vv).TransformationValues(2:3:24)-vmr.Trf(vv).TransformationValues(2));
        % swap superior and inferior coordinates
        tmp=vmr.Trf(vv).TransformationValues(17);
        vmr.Trf(vv).TransformationValues(17)=vmr.Trf(vv).TransformationValues(17-3);
        vmr.Trf(vv).TransformationValues(17-3)=tmp;
        clear tmp;
      else
        error('the transformation matrix in %s%s is not compatible wiht this function.',vmrfname,vmrext);
      end
    end
  end

  if XYZ_flg(3) % along left-right
    try
      vmr.VMRData=flipdim(vmr.VMRData,3);     %#ok
      vmr.VMRData16=flipdim(vmr.VMRData16,3); %#ok
    catch
      vmr.VMRData=flip(vmr.VMRData,3);
      vmr.VMRData16=flip(vmr.VMRData16,3);
    end

    % flip slice position
    vmr.OffsetZ=-1*vmr.OffsetZ;
    vmr.Slice1CenterZ=-1*vmr.Slice1CenterZ;
    vmr.SliceNCenterZ=-1*vmr.SliceNCenterZ;
    if ~overwrite_flg, savefname=[savefname,'Z']; end %#ok

    % flip the transformation matrix
    for vv=1:1:size(vmr.Trf,2)
      if vmr.Trf(vv).NrOfSpatialTransformationValues==16 && ( size(vmr.Trf(vv).TransformationValues,1)==4 && size(vmr.Trf(vv).TransformationValues,2)==4 )
        vmr.Trf(vv).TransformationValues=[1,1,-1,1; 1,1,-1,1; -1,-1,1,-1; 1,1,1,1].*vmr.Trf(vv).TransformationValues;
      elseif vmr.Trf(vv).NrOfSpatialTransformationValues==24 && size(vmr.Trf(vv).TransformationValues,1)==24 % Talairach transformation
        vmr.Trf(vv).TransformationValues(3:3:24)=...
            vmr.Trf(vv).TransformationValues(3)-(vmr.Trf(vv).TransformationValues(3:3:24)-vmr.Trf(vv).TransformationValues(3));
        % swap left and right coordinates
        tmp=vmr.Trf(vv).TransformationValues(24);
        vmr.Trf(vv).TransformationValues(24)=vmr.Trf(vv).TransformationValues(24-3);
        vmr.Trf(vv).TransformationValues(24-3)=tmp;
        clear tmp;
      else
        error('the transformation matrix in %s%s is not compatible wiht this function.',vmrfname,vmrext);
      end
    end
  end
  fprintf('done.\n');

  % setup for a special flipping case -- left/right flip
  if XYZ_flg(1)==0 && XYZ_flg(2)==0 && XYZ_flg(3)==1
    if ~overwrite_flg, savefname=[vmrfname,'_flipLR']; end
  end

  fprintf('saving    : %s%s...',savefname,vmrext);
  vmr.SaveAs(fullfile(vmrpath,[savefname,vmrext]));
  fprintf('done.\n');

  vmr.ClearObject(); clear vmr;
end

return

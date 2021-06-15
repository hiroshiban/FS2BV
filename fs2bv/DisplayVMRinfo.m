function DisplayVMRinfo(VMR_dir,prefix_vmr,save_flg)

% Displays VMR information.
% function DisplayVMRinfo(VMR_dir,:prefix_vmr,:save_flg)
% (: is optional)
%
% This function displays VMR information on the MATLAB command window
%
% [input]
% VMR_dir    : Target directory that contains VMR files
%              e.g. '\CD\3d'
%              Target directory should be specified as such
%              the current directory where this function is
%              called is the origin.
% prefix_vmr : (optional) string to determine the target from
%              multiple files, e.g. 'CD'
% save_flg   : (optional) whether save the displayed info to text file
%              [0|1], 0 by default
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
% Created    : "2012-12-07 17:05:00 banh"
% Last Update: "2018-09-03 12:30:46 ban"

% check input variables
if nargin<1 || isempty(VMR_dir), help(mfilename()); return; end
if nargin<2 || isempty(prefix_vmr), prefix_vmr=''; end
if nargin<3 || isempty(save_flg), save_flg=0; end

if ~exist(fullfile(pwd,VMR_dir),'dir')
  error('can not find VMR_dir. check input variable.');
end

% get VMR file(s)
vmrfiles=GetFiles(fullfile(pwd,VMR_dir),'*.vmr',prefix_vmr);
if isempty(vmrfiles)
  error('VMR file not found. check input variable.');
end

% processing

if save_flg
  fid=fopen(fullfile(pwd,'VMRinfo.txt'),'w');
  if fid==-1, error('can not open VMRinfo.txt to write.'); end
end

fprintf('\nDisplaying VMR file information\n\n');
fprintf('Target directory: %s\n\n',fullfile(pwd,VMR_dir));
if save_flg
  fprintf(fid,'Target directory: %s\n\n',fullfile(pwd,VMR_dir));
end

vmrcounter=0;
for ii=1:1:length(vmrfiles)
  [vmrpath,vmrname,vmrext]=fileparts(vmrfiles{ii});

  if ~strcmpi(vmrext,'.vmr') % to avoid some file like *.vmr_3DMC_verbose.log
    continue
  else
    vmrcounter=vmrcounter+1;
  end

  fprintf('VMR #%02d: %s%s\n',vmrcounter,vmrname,vmrext);
  if save_flg, fprintf(fid,'VMR #%02d: %s%s\n',vmrcounter,vmrname,vmrext); end

  vmr=BVQXfile(vmrfiles{ii});

  fprintf('\n');
  fprintf('     FOV [row,col] = [%d, %d] (mm)\n',vmr.FoVRows,vmr.FoVCols);
  fprintf('     data dimension [x,y,z] = [%d, %d, %d] (voxels)\n',vmr.DimX,vmr.DimY,vmr.DimZ);
  fprintf('     offset [x,y,z] = [%d, %d, %d] (voxels)\n',vmr.OffsetX,vmr.OffsetY,vmr.OffsetZ);
  fprintf('     voxels [x,y,z] = [%.2f, %.2f, %.2f] (mm)\n',vmr.VoxResX,vmr.VoxResY,vmr.VoxResZ);
  fprintf('\n');
  fprintf('     first slice image center [x,y,z] = [%.4f, %.4f, %.4f]\n',...
          vmr.Slice1CenterX,vmr.Slice1CenterY,vmr.Slice1CenterZ);
  fprintf('     last slice image center [x,y,z] = [%.4f, %.4f, %.4f]\n',...
          vmr.SliceNCenterX,vmr.SliceNCenterY,vmr.SliceNCenterZ);
  fprintf('     raw direction [x,y,z] = [%.4f, %.4f, %.4f]\n',...
          vmr.RowDirX,vmr.RowDirY,vmr.RowDirZ);
  fprintf('     col direction [x,y,z] = [%.4f, %.4f, %.4f]\n',...
          vmr.ColDirX,vmr.ColDirY,vmr.ColDirZ);
  fprintf('\n');

  if save_flg
    fprintf(fid,'\n');
    fprintf(fid,'     FOV [row,col] = [%d, %d] (mm)\n',vmr.FoVRows,vmr.FoVCols);
    fprintf(fid,'     data dimension [x,y,z] = [%d, %d, %d] (voxels)\n',vmr.DimX,vmr.DimY,vmr.DimZ);
    fprintf(fid,'     offset [x,y,z] = [%d, %d, %d] (voxels)\n',vmr.OffsetX,vmr.OffsetY,vmr.OffsetZ);
    fprintf(fid,'     voxels [x,y,z] = [%.2f, %.2f, %.2f] (mm)\n',vmr.VoxResX,vmr.VoxResY,vmr.VoxResZ);
    fprintf(fid,'\n');
    fprintf(fid,'     first slice image center [x,y,z] = [%.4f, %.4f, %.4f]\n',...
            vmr.Slice1CenterX,vmr.Slice1CenterY,vmr.Slice1CenterZ);
    fprintf(fid,'     last slice image center [x,y,z] = [%.4f, %.4f, %.4f]\n',...
            vmr.SliceNCenterX,vmr.SliceNCenterY,vmr.SliceNCenterZ);
    fprintf(fid,'     raw direction [x,y,z] = [%.4f, %.4f, %.4f]\n',...
            vmr.RowDirX,vmr.RowDirY,vmr.RowDirZ);
    fprintf(fid,'     col direction [x,y,z] = [%.4f, %.4f, %.4f]\n',...
            vmr.ColDirX,vmr.ColDirY,vmr.ColDirZ);
    fprintf(fid,'\n');
  end

  vmr.ClearObject(); clear vmr;
end

if save_flg, fclose(fid); end

return

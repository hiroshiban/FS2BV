function LimitVOIvoxelsByVMPthreshold(VOI_dir,vmp_file,mapnum,ignore_sign_flg,thres,unique_flg,prefix_voi,overwrite_flg)

% Limits the number of voxels in each VOI based on the specified VMP response map (e.g. t-values).
% function LimitVOIvoxelsByVMP(VOI_dir,vmp_file,:mapnum,:ignore_sign_flg,:thres,:unique_flg,:prefix_voi,:overwrite_flg)
% (: is optional)
%
% Limits the number of voxels in VOI file(s) by the values in vmp file
%
% [input]
% VOI_dir    : Target directory that contains VOI files
%              e.g. '\CD\zk08_382\voi_files'
%              Target directory should be specified as such
%              the current directory where this function is
%              called is the origin.
% vmp_file   : vmp file to be used for thresholding voxels
%              specify with relative path from the current directory
%              e.g. vmp_file='./zk11_052/ROI_vtc/vmp/zk11_052.vmp';
% mapnum     : (optional) map id in the vmp_file, 1 by default
% ignore_sign_flg : (optional) whether ignoring the sign of the VMP voxel values.
%              true (ignore the sign and use the absolute values) or false. false by default.
% thres      : (optional) a threshold, voxels with the values above this
%              will be survived. 2.0 by default
% unique_flg : whether using only unique voxel VMP values or not, 1 by default
%              if 0, a same voxel value may be counted more than once
%              (in same cases, this will be valid)
% prefix_voi : (optional) string to determine the target VOI
%              from multiple files, e.g. '*combined_all'
% overwrite_flg: (optional) whether overwrite the original file with the new one.
%                if set to 1, all the files will be overwritten by new ones.
%                if set to 0, the new data will be saved with '_limit_thr{val}'
%                prefix in the same directory with the input VOIs. 0 by default.
%
% [output]
% no output variable,
% the thresholded VOIs overwrite the existing VOIs (when overwrite_flg is 1)
% or will be saved as (voi_name)_limit{nvox}.voi (when overwrite_flg is 0).
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
% Created    : "2012-03-15 12:10:33 banh"
% Last Update: "2018-09-03 11:33:28 ban"

% check input variable
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(mapnum), mapnum=1; end
if nargin<4 || isempty(ignore_sign_flg), ignore_sign_flg=false; end
if nargin<5 || isempty(thres), thres=2.0; end
if nargin<6 || isempty(unique_flg), unique_flg=1; end
if nargin<7 || isempty(prefix_voi), prefix_voi=''; end
if nargin<8 || isempty(overwrite_flg), overwrite_flg=0; end

% get VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if isempty(voifiles)
  error('no voi files found. check input variable');
end

% get the target VMP file
[rpath,vmpname,vmpext]=fileparts(vmp_file);
vmpfile=GetFiles(fullfile(pwd,rpath),['*',filesep(),vmpname,vmpext],'',1,1);
if isempty(vmpfile)
  error('target vmp file not found. check input variable.');
end

% load VMP file and store the voxel information
vmp=BVQXfile(vmpfile{1});
if mapnum>vmp.NrOfMaps
  error('the target vmp only has %d maps (<%d), check input variable.',vmp.NrOfMaps,mapnum);
end
vox_resolution=vmp.Resolution;
Xs=vmp.XStart;
Ys=vmp.YStart;
Zs=vmp.ZStart;
%Xe=vmp.XEnd;
%Ye=vmp.YEnd;
%Ze=vmp.ZEnd;
vmpData=vmp.Map(mapnum).VMPData;
vox_dims=size(vmpData);

% display message
fprintf('target VOI dir : %s\n',fullfile(pwd,VOI_dir));
fprintf('target VMP file: %s%s\n',vmpname,vmpext);
fprintf('\n');

% processing
for ii=1:1:length(voifiles)

  % display message
  [path,voifile,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s\n',voifile,voiext);

  % load VOI
  voi=BVQXfile(voifiles{ii});
  org_VMRdims=voi.OriginalVMRFramingCubeDim;

  % extract VOI voxel coordinates
  fprintf('extracting VOI voxels...');
  voi_idx=cell(voi.NrOfVOIs,1);
  for vv=1:1:voi.NrOfVOIs
    voi_voxels=voi2col(voi.VOI(vv).NrOfVoxels,voi.VOI(vv).Voxels,vox_resolution,'off',[Xs,Ys,Zs],vox_dims,org_VMRdims,voi.ReferenceSpace);
    voi_idx{vv}=sub2ind(vox_dims,voi_voxels(:,1),voi_voxels(:,2),voi_voxels(:,3));
  end
  fprintf('done.\n');

  % thresholding
  fprintf('thresholding VOI voxels...');
  voxels=cell(voi.NrOfVOIs,1);
  for vv=1:1:voi.NrOfVOIs
    mapvals=vmpData(voi_idx{vv});
    if ignore_sign_flg, mapvals=abs(mapvals); end

    % if you want to use only unique voxel t-value
    if unique_flg
      [mapvals,m]=unique(mapvals);
      unique_idx=voi_idx{vv}(m);

      [dummy,order_idx]=sort(mapvals,'descend'); %#ok
      cvox_idx=unique_idx(order_idx);
    else
      [dummy,order_idx]=sort(mapvals,'descend'); %#ok
      cvox_idx=voi_idx{vv}(order_idx);
    end

    % thresholding VMP voxels
    if ignore_sign_flg
      cvox_idx=cvox_idx(abs(vmpData(cvox_idx))>thres);
    else
      cvox_idx=cvox_idx(vmpData(cvox_idx)>thres);
    end

    % store (x,y,z) coordinates of voxels in the whole Tal space
    if isempty(cvox_idx)
      voxels{vv}(1,1:3)=NaN;
    else
      [voxels{vv}(:,1),voxels{vv}(:,2),voxels{vv}(:,3)]=ind2sub(vox_dims,cvox_idx);
    end
  end
  fprintf('done.\n');

  % updating thresholded voxels to VOI structures
  fprintf('updating thresholded voxel information...');

  if vox_resolution==1 || vox_resolution==2
    offset=[1,1,1];
  elseif vox_resolution==3
    offset=[2,2,2];
  else
    error('VMP''s voxel resolution should be one of 1,2,3.');
  end

  for vv=1:1:voi.NrOfVOIs
    voi.VOI(vv).NrOfVoxels=size(voxels{vv},1);
    tmp_vox=ones(size(voxels{vv},1),1)*ceil( 128-([Xs,Ys,Zs]-offset) )-voxels{vv}.*vox_resolution;
    voi.VOI(vv).Voxels=tmp_vox(:,[3,1,2]);
  end
  fprintf('done.\n');

  fprintf('saving the new voi...');
  if overwrite_flg
    voi.SaveAs(fullfile(path,[voifile,voiext]));
  else
    voi.SaveAs(fullfile(path,[voifile,sprintf('_limit_thr%.2f',thres),voiext]));
  end
  fprintf('done.\n');

  % clean up
  voi.ClearObject(); clear voi;

end % for ii=1:1:length(voifiles)

% clean up
vmp.ClearObject(); clear vmp;

return

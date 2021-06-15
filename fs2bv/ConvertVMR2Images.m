function ConvertVMR2Images(VMR_dir,oimg_name,oimg_format,prefix_vmr,orientation)

% Convertes VMR anatomy data to multiple 2D image files (bmp,gif,jpg,png,tiff).
% function ConvertVMR2Images(VMR_dir,oimg_name,:oimg_format,:prefix_vmr,:orientation)
% (: is optional)
%
% [about]
% This simple function converts VMR files to multiple image files
%
% [example]
% >> ConvertVMR2Images('/AA/zk08_382/','AA_3D','bmp','zk')
%
% [input]
% VMR_dir      : Target directory that contains VMR files
%                e.g. '\CD\zk08_382'
%                Target directory should be specified as such
%                the current directory where this function is
%                called is the origin.
%
% oimg_name    : output file name, string
%                e.g. if oimg_name='3D' & oimg_format='bmp', then
%                     3D_001.bmp, 3D_002.bmp, ... and so on are stored
%                     in VMR_dir/Imgs
%
% oimg_format  : (optional) output file format
%                one of 'bmp', 'gif', 'jpg', 'png', or 'tiff'
% prefix_vmr   : (optional) string to determine the target VMR
%                from multiple files, e.g. 'CD'
% orientation  : (optional) orientation of the converted image(s).
%                1: saggittal
%                2: coronal
%                3: axial
%                NOTE: the original VMRs should be in the default
%                SAG coordinate of BrainVoyager QX
%
% [output]
% No output variable, converted images are strored in VMR_dir/Imgs
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
% Created    : "2010-08-18 11:50:53 banh"
% Last Update: "2018-09-03 10:53:25 ban"

% check the input variables
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(oimg_format), oimg_format='bmp'; end
if nargin<4 || isempty(prefix_vmr), prefix_vmr=''; end
if nargin<5 || isempty(orientation), orientation=1; end

% get vmr files from the directories in the VMR_dir
vmrfiles=GetFiles(fullfile(pwd,VMR_dir),'*.vmr',prefix_vmr);

% display message
fprintf('Target : %s\n',fullfile(pwd,VMR_dir));

for ii=1:1:length(vmrfiles)

  % display message
  [path,file,ext]=fileparts(vmrfiles{ii});
  fprintf('processing : %s%s\n',file,ext);

  % create Imgs directory
  if ~exist(fullfile(path,'Imgs'),'dir'), mkdir(fullfile(path,'Imgs')); end

  % read vtc file & create its object
  vmr=BVQXfile(vmrfiles{ii});

  % rotate vmr anatomical data
  % this rotation is compatible with BVQX coordinate,
  % but please be careful and change this line if you
  % want to get images from different orientation
  if orientation==1
    vmr.VMRData=permute(vmr.VMRData,[2 1 3]); % saggittal
  elseif orientation==2
    vmr.VMRData=permute(vmr.VMRData,[2 3 1]); % coronal
  elseif orientation==3
    vmr.VMRData=permute(vmr.VMRData,[1 3 2]); % axial
  else
    error('orientation should be one of 1-3. check input variables.');
  end

  for zz=1:1:size(vmr.VMRData,3)
    imname=sprintf([oimg_name,'_%03d'],zz);
    imwrite(vmr.VMRData(:,:,zz),[fullfile(path,'Imgs'),filesep(),imname,'.',oimg_format],oimg_format);
  end

  % clear BVQX object
  vmr.ClearObject(); clear vmr;

end % for ii=1:1:length(vtcfiles)

return

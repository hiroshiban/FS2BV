function OverlayVOI2VMR(VMR_dir,VOI_dir,output_format,orientation,voi_to_use,alphaval,colors,prefix_vmr,prefix_voi,border_flg)

% Overlays VOIs on VMP anatomy file and saves the overlayed data as 2D images.
% function OverlayVOI2VMR(VMR_dir,VOI_dir,:output_format,:orientation,:voi_to_use,:alphaval,:colors,:prefix_vmr,:prefix_voi,:border_flg)
% (: is optional)
%
% Overlays VOIs on VMR and saves generated images as output_format.
%
% [example]
% >> VMR_dir='HB/3d'; VOI_dir='HB/zk10_148/voi_files/
% >> output_format='bmp';
% >> orientation='SAG';
% >> voi_to_use={'V1','V2','V3','V3A'};
% >> alphaval=0.4;
% >> OverlayVOI2VMR(VMR_dir,VOI_dir,output_format,orientation,voi_to_use,alphaval,[]);
%
% [input]
% VMR_dir   : Target directory that contains VMR file
%             e.g. '\HB\zk10_148\3d'
%             Target directory should be specified as such
%             the current directory where this function is
%             called is the origin.
% VOI_dir   : Target directory that contains VMR file
%             e.g. '\HB\zk10_148\3d'
%             Target directory should be specified as such
%             the current directory where this function is
%             called is the origin.
% output_format : (optional) output image format, one of
%                 'BMP','GIF','HDF','JPEG','PBM','PCX','PNG','PPM','RAS','TIFF','XWD'.
%                 'BMP' by default
% orientation   : output image orientation, one of 'SAG', 'COR', 'TRA'. 'SAG' by default
% voi_to_use    : names of VOIs to be overlayed, cell structure, e.g. {'V1','V2','V3A'}
%                 If voi_to_use=[] (empty), all the VOIs in VOI files will be overlayed.
% alphaval      : alpha of overlay, [0.0-1.0]. 0.4 by default
% colors        : (optional) VOI colors you want to use
%                 n x 3 matrix, e.g. colors=[1,0,0;0,1,0;...]; (num_VOI x RGB(0.0-1.0))
%                 if set to 0, original colors in VOI files will be used.
%                 if set to empty, rainbow colors will be used. colors=[]; by default.
%                 if size(colors,1)<length(voi_to_use), the final RGB color will be used
%                 for coloring the rest of VOIs
% prefix_vmr    : (optoinal) string to determine the target from
%                 multiple files, e.g. '*_final_TAL'
%                 '*_final_TAL' by default
% prefix_voi    : (optoinal) string to determine the target from
%                 multiple files, e.g. 'all_ROIs_combined*', empty by default
% border_flg    : (optional) if 1, only the border of each of VOIs is painted.
%                 0 by default.
%
% [output]
% no output variable
% generated images will be saved in ~/vmr_voi_imgs
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
% Created    : "2012-11-23 01:32:57 banh"
% Last Update: "2021-04-28 18:27:24 ban"

% check input variables
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(output_format), output_format='BMP'; end
if nargin<4 || isempty(orientation), orientation='SAG'; end
if nargin<5 || isempty(voi_to_use), voi_to_use={}; end
if nargin<6 || isempty(alphaval), alphaval=0.4; end
% "colors" variable will be set later
if nargin<8 || isempty(prefix_vmr), prefix_vmr='*_final_TAL'; end
if nargin<9 || isempty(prefix_voi), prefix_voi=''; end
if nargin<10 || isempty(border_flg), border_flg=0; end

if ~iscell(voi_to_use), voi_to_use={voi_to_use}; end

% get VMR file
if ~exist(fullfile(pwd,VMR_dir),'dir')
  error('VMR_dir not found. check input variable');
end
vmrfiles=GetFiles(fullfile(pwd,VMR_dir),'*.vmr',prefix_vmr);
if length(vmrfiles)>1
  fprintf('multiple VMR files are found..using the first one...\n');
  vmrfiles={vmrfiles{1}};
end
[vmrpath,vmrname,vmrext]=fileparts(vmrfiles{1});
fprintf('TARGET DIR : %s\n',vmrpath);
fprintf('       VMR : %s%s\n',vmrname,vmrext);

% get VOI file(s)
if ~exist(fullfile(pwd,VOI_dir),'dir')
  error('VOI_dir not found. check input variable');
end
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if length(voifiles)>1
  CombineVOIs(VOI_dir,'tmp_combined.voi',prefix_voi);
  voifile=fullfile(pwd,VOI_dir,'tmp_combined.voi');
else
  voifile=voifiles{1};
end
[voipath,voiname,voiext]=fileparts(voifile);
fprintf('TARGET DIR : %s\n',voipath);
fprintf('       VOI : %s%s\n',voiname,voiext);

% processing VMR & VOI
fprintf('loading VMR...');
vmr=BVQXfile(vmrfiles{1});
VMRData=vmr.VMRData;
vmr.ClearObject(); clear vmr;
fprintf('done.\n');

fprintf('loading target VOIs...');
voi=BVQXfile(voifile);
idx=0;
if ~isempty(voi_to_use)
  VOIs=cell(length(voi_to_use),1);
  for vv1=1:1:length(voi_to_use)
    preidx=idx;
    for vv2=1:1:voi.NrOfVOIs
      if strcmp(voi.VOI(vv2).Name,voi_to_use{vv1})
        idx=idx+1;
        VOIs{idx}=voi.VOI(vv2);
      end
    end
    if idx==preidx
      fprintf('VOI: %s not found. skipping...\n',voi_to_use{vv1});
    end
  end
else
  VOIs=cell(length(voi.NrOfVOIs),1);
  for vv2=1:1:voi.NrOfVOIs
    VOIs{vv2}=voi.VOI(vv2);
  end
end
fprintf('done.\n');

%% overlaying VOIs on VMR

% converting 1D-grayscale image to 3D-RGB image
VMR=reshape(repmat(VMRData(:),[1,3]),[size(VMRData),3]);

% set colors
if nargin<7 || isempty(colors) % use rainbow colors
  colors=rainbow(length(VOIs));
elseif colors==0 % use original colors described in VOI file
  colors=zeros(length(VOIs),3);
  for vv=1:1:length(VOIs), colors(vv,:)=VOIs{vv}.Color; end
end
if max(colors(:))<=1.0, colors=ceil(colors.*255); end

% overlapping VOIs on VMR
fprintf('overlaying VOIs on VMR images...');
for vv=1:1:length(VOIs)
  vox_ind=-VOIs{vv}.Voxels+128;

  if border_flg
    % check whether the target voxel is on the VOI border or not based on the neighboring 26 values of the target voxel
    % -- if neighbors contain 1, the target voxel is regarded to be on the border
    bvox_ids=zeros(size(vox_ind,1),1);
    bbox=ones(size(VMRData)); % a box to store the outside regions of the VOI
    bbox(sub2ind(size(VMRData),vox_ind(:,2),vox_ind(:,3),vox_ind(:,1)))=0; % 0 in the target VOI region
    bbox_updated=bbox;
    for mm=1:1:2 % repeat voxel selection to make the line thicker
      for ii=1:1:size(vox_ind,1)
        neighbors=bbox(vox_ind(ii,2)-1:vox_ind(ii,2)+1,vox_ind(ii,3)-1:vox_ind(ii,3)+1,vox_ind(ii,1)-1:vox_ind(ii,1)+1);
        neighbors=neighbors(:);
        neighbors(14)=[];

        if sum(neighbors)>8 %if find(neighbors==1,1)
          bvox_ids(ii)=1;
          bbox_updated(vox_ind(ii,2),vox_ind(ii,3),vox_ind(ii,1))=1;
        end
      end
      bbox=bbox_updated;
    end

    % update vox_ind so that it only contains the border voxels
    vox_ind=vox_ind(logical(bvox_ids),:);
  end

  for ii=1:1:size(vox_ind,1)
    for cc=1:1:3 % RGB
      VMR(vox_ind(ii,2),vox_ind(ii,3),vox_ind(ii,1),cc)=...
          (1-alphaval)*VMR(vox_ind(ii,2),vox_ind(ii,3),vox_ind(ii,1),cc)+...
          alphaval*colors(min(vv,size(colors,1)),cc);
    end
  end
end
fprintf('done.\n');
voi.ClearObject(); clear voi;

% set image orientation
if strcmpi(orientation,'SAG')
  VMR=permute(VMR,[2,1,3,4]); % sagittal
elseif strcmpi(orientation,'COR')
  VMR=permute(VMR,[2,3,1,4]); % colonal
elseif strcmpi(orientation,'TRA')
  VMR=permute(VMR,[1,3,2,4]); % axial
else
  error('orientation should be one of ''SAG'',''COR'',''TRA''. check input variable.');
end

% saving data
fprintf('saving data as %s format...',lower(output_format));
save_dir=fullfile(pwd,'vmr_voi_imgs');
if ~exist(save_dir,'dir'), mkdir(save_dir); end

f1=figure('Name','VMR images overlayed with VOIs','NumberTitle','off','MenuBar','none','ToolBar','none');
for ii=1:1:size(VMR,3), imshow(squeeze(VMR(:,:,ii,:))); drawnow; end
for ii=1:1:size(VMR,3)
  fname=sprintf('img_%s_%03d.%s',orientation,ii,lower(output_format));
  imwrite(squeeze(VMR(:,:,ii,:)),fullfile(save_dir,fname),lower(output_format));
end
close(f1);
fprintf('done.\n');

return


%% subfunction

function map = rainbow(m)
%
%   function map = rainbow(m)
%
%   RAINBOW(M) returns an M-by-3 matrix containing an RAINBOW colormap.
%   RAINBOW, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(rainbow)
%
%   See also GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM, JET,
%   COLORMAP, RGBPLOT, HSV2RGB, RGB2HSV.
%
% March 98H.Yamamoto
% Modified by Hiroshi Ban, March 08 2011

if nargin<1, m = size(get(gcf,'colormap'),1); end
% h = (0:m-1)'/max(m,1);
h = (m-1:-1:0)'/max(m,1);
if isempty(h)
  map = [];
else
  map = hsv2rgb([h ones(m,2)]);
end

return

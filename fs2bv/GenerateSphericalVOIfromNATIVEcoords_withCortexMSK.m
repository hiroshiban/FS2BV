function voi=GenerateSphericalVOIfromNATIVEcoords_withCortexMSK(SRF_dir,center_XYZ,radius_mm,voinames,voicolors,save_name,prefix_srf,srf_params,vmr_params,VMRdim)

% Generates spherical VOIs in the NATIVE space, and masks them with a cortex mask.
% function GenerateSphericalVOIfromNATIVEcoords_withCortexMSK(SRF_dir,center_XYZ,radius_mm,voinames,:voicolors,:save_name,:prefix_srf,:srf_params,:vmr_params,:VMRdim)
% (: is optional)
%
% This function generates spheric VOI(s) as the center='center_XYZ' with radius='radius_mm'
% And it masks the generated spheric VOI by cortex graymatter mask generated from *RECOSM.SRF_{LH|RH}.srf
%
% [example]
% >> center_XYZ=[-35.3,-80.8,5.0;31.7,-84.7,2.5];
% >> voinames={'KO_left_tal','KO_right_tal'};
% >> save_name='KO_talairach.voi';
% >> SRF_dir='/HB/3d/';
% >> voi=GenerateSphericalVOIfromNATIVEcoords_withCortexMSK(SRF_dir,center_XYZ,5,voinames,[],save_name,'*RECOSM');
%
% [input]
% SRF_dir         : target directory that contains RECOSM SRF file (*RECOSM_{LH|RH}.srf)
%                   e.g. '/HB/3D/zk09_034/
%                   the location should be specified as such
%                   the current directory where this function is
%                   called is the origin.
% center_XYZ  : talairach coordinate, num_voi x 3 matrix, [X1,Y1,Z1; X2,Y2,Z2,...]
% radius_mm       : radius of the ROI in mm, [val]
% voinames        : names of VOIs, cell structure, {'ROI_1','ROI_2',...}
% voicolors       : (optional) colors of VOIs, cell structure, {[red1,green1,blue1],[red2,...]}
% save_name       : (optional) if not empty, the generated voi structure will be saved
%                   with this name in the current directory, 'my_vois.voi'
% prefix_srf      : (optoinal) string to determine the target srf
%                   from multiple files, e.g. 'HB'
% srf_params      : (optional) parameters used for srf.Combine (combine left/rigth srfs)
%   parameters are as below
%     srf_params
%       .type       one of
%                   'backtoback' - rotate one mesh 180 in XY plane
%                   'custom'     - build custom scene, see below
%                   'gapped'     - join contents with a gap
%                   'wholebrain' - simply join contents (default)
%       .color1     1x4 double for .Color1 field in SRF options
%       .color2     1x4 double for .Color2 field in SRF options
%       .filename   store combined under new filename
%       .gap        1x1 double, mm to insert between two meshes
%                   (applied for all types accordingly, defaults:
%                   backtoback: 25, gapped: 100, outandin: 20,
%                   outintb: 25, patched: 25, spm2: 25, wholebrain: 0)
%       .linkedmtc  filename of linked MTC, set to empty if not given
%       .mtc1       BVQXfile MTC object for the first SRF
%       .mtc2       BVQXfile MTC object for the second SRF
%       .smp1       BVQXfile SMP object for the first SRF
%       .smp2       BVQXfile SMP object for the second SRF
%       .ssm1       BVQXfile SSM object for the first SRF
%       .ssm2       BVQXfile SSM object for the second SRF
%       .transform  1x2 cell array with 4x4 double transformation matrices,
%                   needed for custom scenaries
% vmr_params      : (optional) parameters used for srf.BackToVMR
%   parameters are as below
%     vmr_params
%        .fillmode  either of {'nearest'} or 'linear'
%        .from      from V + from * N (default: -0.75)
%        .step      step (default: 0.5)
%        .to        to V + to * N (default: 0.75)
% VMRdim          : (optional) VMR volume dimensions in which the generated VOIs
%                   to be projected. 256 by default (1mm-cubic VMR), use 256 for 0.5mm VMR.
%
% [output]
% voi             : generated VOI structure, BVQX voi data
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
% Created    : "2020-03-30 16:44:32 ban"
% Last Update: "2020-03-30 17:20:38 ban"

%% check the input variables
if nargin<4, help(mfilename()); return; end

if ~exist(fullfile(pwd,SRF_dir),'dir')
  error('SRF_dir is not found. check input variable.');
end

if nargin<5 || isempty(voicolors)
  voicolors=cell(size(center_XYZ,1),1);
  tcolors=rainbow(size(center_XYZ,1));
  for ii=1:1:size(center_XYZ,1)
    voicolors{ii}=tcolors(ii,:);
  end
end

if nargin<6, save_name=''; end
if ~isempty(save_name) && ~strcmp(save_name(end-3:end),'.voi')
  save_name=strcat(save_name,'.voi');
end

if nargin<8 || isempty(prefix_srf), prefix_srf=''; end

% check parameters for combining left/right srfs
if nargin<9 || isempty(srf_params)
  srf_params.type='wholebrain';
  srf_params.gap=0;
else
  if ~isstructmember(srf_params,'type')
    srf_params.type='wholebrain';
  end
  if ~isstructmember(srf_params,'gap')
    srf_params.gap=0;
  end
end

% check parameters for converting from srf to vmr(msk)
if nargin<10 || isempty(vmr_params)
  vmr_params.fillmode='nearest';
  vmr_params.from=-1;
  vmr_params.step=0.5;
  vmr_params.to=3;
else
  if ~isstructmember(vmr_params,'fillmode')
    vmr_params.fillmode='nearest';
  end
  if ~isstructmember(vmr_params,'from')
    vmr_params.from=-1;
  end
  if ~isstructmember(vmr_params,'step')
    vmr_params.step=0.5;
  end
  if ~isstructmember(vmr_params,'to')
    vmr_params.to=3;
  end
end

if nargin<11 || isempty(VMRdim), VMRdim=256; end

if size(center_XYZ,2)~=3
  error('center_XYZ should be [npoints x 3] matrix. check input variable');
end

if size(center_XYZ,1)~=length(voinames)
  error('the numbers of center_XYZ and voinames mismatch. check input variable');
end

if ~iscell(voinames), voinames={voinames}; end
if ~iscell(voicolors), voicolors={voicolors}; end

scaling=VMRdim/256;

%%% generate cortex graymatter mask from *RECOSM.srf files
fprintf('Generating Graymatter mask...\n');

% constant variable, mask_resolution should be always 1 to generate VOIs
mask_resolution=1;

% get srf files from the directories in the SRF_dir
srffiles=GetFiles(fullfile(pwd,SRF_dir),'*.srf',prefix_srf);

% get recosm files from the directories in the RECOSM_dir
if length(srffiles)>3, error('more than 2 SRF files are found! Specify prefix_srf and try again'); end

% display message
fprintf('Target SRF   : %s\n',fullfile(pwd,SRF_dir));

% set left or right prefix for RECOSM files
leftcounter=0; rightcounter=0;
for ii=1:1:length(srffiles)
  [tpath,srfname]=fileparts(srffiles{ii}); %#ok
  if strfind(srfname,'LH')
    leftcounter=leftcounter+1;
    srfleft{leftcounter}=srffiles{ii};  %#ok
  elseif strfind(srfname,'RH')
    rightcounter=rightcounter+1;
    srfright{rightcounter}=srffiles{ii}; %#ok
  end
end

% check whether existing multiple files in one hemisphere
if length(srfleft)>=2 || length(srfright)>=2
  error('Multiple reference srf files in one hemisphere. Specify prefix_srf and try again');
end

% get path of left srf
[path,srfname,ext]=fileparts(srfleft{1}); %#ok

% display message
fprintf('processing : %s%s\n',srfname,ext);

% read srf file & create its object
srf_LH=BVQXfile(srfleft{1});
srf_RH=BVQXfile(srfright{1});

% check whether reference file is attached to the original srfs
if ( isempty(strfind(srfleft{1},'RECO')) && isempty(strfind(srfleft{1},'RECOSM')) )
  if isempty(srf_LH.AutoLinkedMTC)
    error('reference srf file is not attached to the left srf! check AutoLinkedMTC');
  else
    new_LH=srf_LH.AutoLinkedMTC;
    srf_LH.ClearObject(); clear srf_LH;
    srf_LH=BVQXfile(new_LH);
  end
end

if ( isempty(strfind(srfright{1},'RECO')) && isempty(strfind(srfright{1},'RECOSM')) )
  if isempty(srf_RH.AutoLinkedMTC)
    error('reference srf file is not attached to the right srf! check AutoLinkedMTC');
  else
    new_RH=srf_RH.AutoLinkedMTC;
    srf_RH.ClearObject(); clear srf_RH;
    srf_RH=BVQXfile(new_RH);
  end
end

% attach the right srf to the left srf.
[path,srfname,ext]=fileparts(srfright{1}); %#ok
fprintf('             <-- %s%s\n',srfname,ext);

srf_LR=srf_LH.Combine(srf_RH,srf_params);

% generate vmr files
vmr=srf_LR.BackToVMR(vmr_params);

% convert vmr to msk
mask=zeros(VMRdim,VMRdim,VMRdim);
mask=mask(vmr.VertexVMRData~=0)=1;
mask=permute(mask,[3,1,2]);

% clean up
clear vmr;
srf_LH.ClearObject(); clear srf_LH;
srf_RH.ClearObject(); clear srf_RH;

% get mask id
idxgray=find(mask~=0);

fprintf('done.\n');

%% generate VOIs

% display message
fprintf('\n');
if size(center_XYZ,1)==1
  fprintf('Generating spheric VOI (radius: %02dmm) in the native coordinates\n',radius_mm);
else
  fprintf('Generating spheric VOIs (radius: %02dmm) in the native coordinates\n',radius_mm);
end

% convert unit of radius
radius_mm=(scaling*radius_mm)^2;

% convert tal to BVQX's voi coords
fprintf('Initializing voi full coordinates...');

% Initially we are presented with a list of voi coordinates
% in the native space. These coordinates range from -128 to 128.
center_voi_XYZ=scaling.*center_XYZ;

% generate 1mm full-space VOI coordinate
% here, 'repmat' more than 1 repetition makes the processing slow
% using 'for' loop partly is faster
voi_full_coords=zeros(VMRdim,VMRdim,VMRdim,3);
coordarray1=repmat((1:VMRdim)'-0.5,1,VMRdim); % -0.5 is required to set the center of voxel
coordarray2=repmat((1:VMRdim)-0.5,VMRdim,1);
for idx=1:1:VMRdim
  voi_full_coords(:,:,idx,1)=coordarray1;
  voi_full_coords(idx,:,:,2)=coordarray1;
  voi_full_coords(idx,:,:,3)=coordarray2;
end
fprintf('done.\n');

% initialize voi object
voi=BVQXfile('new:voi');
voi.FileVersion=4;
%voi.SubjectVOINamingConvention='<VOI>_<SUBJ>';
voi.NrOfVOIs=size(center_voi_XYZ,1);
voi.OriginalVMRResolutionX=256/VMRdim;
voi.OriginalVMRResolutionY=256/VMRdim;
voi.OriginalVMRResolutionZ=256/VMRdim;
voi.OriginalVMRFramingCubeDim=VMRdim;

% setting new VOI details
for ii=1:1:size(center_voi_XYZ,1)

  % display message
  fprintf('processing (NATIVE:% 6.1f,% 6.1f,% 6.1f) as % 14s...',...
          center_XYZ(ii,1),center_XYZ(ii,2),center_XYZ(ii,3),voinames{ii});

  % generate spheric coordinates
  coord1=voi_full_coords(:,:,:,1)-center_voi_XYZ(ii,1);
  coord2=voi_full_coords(:,:,:,2)-center_voi_XYZ(ii,2);
  coord3=voi_full_coords(:,:,:,3)-center_voi_XYZ(ii,3);
  distances=coord1.^2+coord2.^2+coord3.^2;
  idx=find(distances<=radius_mm);

  % mask generated VOI coords by gray-matter mask
  idx=intersect(idx,idxgray);

  % note that BVQX's coords [x,y,z] = [y,z,x] in MATLAB
  [XX,YY,ZZ]=ind2sub([VMRdim,VMRdim,VMRdim],idx);

  % set the coordinates to voi object
  for vv=1:1:numel(XX)
    voi.VOI(ii).Voxels(vv,:)=[VMRdim/2-XX(vv),VMRdim/2-YY(vv),VMRdim/2-ZZ(vv)];
  end

  % set the additional parameters to voi object
  voi.VOI(ii).Name=voinames{ii};
  voi.VOI(ii).Color=voicolors{ii};
  voi.VOI(ii).NrOfVoxels=numel(XX);

  fprintf('done.\n');

end

% save the VOI
if ~isempty(save_name)
  voi.SaveAs(fullfile(pwd,save_name));
end

if ~nargout
  voi.ClearObject(); clear voi;
end

fprintf('completed.\n');

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
  map = floor(255*hsv2rgb([h ones(m,2)]));
end

return

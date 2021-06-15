% script_convert_networks_NII2VOIs.m
%
% [usage]
% 1. please copy networks.nii, networks.txt, and networks.info onto the directory where this script is located
% 2. please run this script on the MATLAB shell.
%
%
% Created    : "2018-08-27 14:41:10 ban"
% Last Update: "2019-12-16 15:14:49 ban"

cv_hbtools_BVQX_setup(1);

% some constants
NII_dir='.';
labels={};

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii='networks';
save_flg=1;

% loading ROI label information
fid=fopen(sprintf('%s.txt',prefix_nii),'r');
roicounter=0;
tmproi=fgetl(fid);
while ischar(tmproi)
  roicounter=roicounter+1;
  ROIs{roicounter}=tmproi; %#ok
  tmproi=fgetl(fid);
end
fclose(fid);
nROIs=length(ROIs);

% set ROI colors
lut=ceil(255.*generate_colormap_rainbow(nROIs));

% set labels for ConvertNiftiRoi2BVvoi_Labels
for ii=1:1:nROIs, labels(ii,:)={ii,ROIs{ii},lut(ii,:)}; end %#ok

% updating the layered nii to indexed nii file
% (the matrix is converted from 91 x 109 x 91 x 32 (32 layers of 0|1) to 91 x 109 x 91).
nii=load_nii(sprintf('%s.nii',prefix_nii));
sz=size(nii.img);
volimg=zeros(sz(1:3));
for vv=1:1:sz(4), volimg(nii.img(:,:,:,vv)==1)=vv; end
nii.img=volimg;
nii.hdr.dime.dim(1)=3;
nii.hdr.dime.dim(5)=1;
nii.hdr.dime.dim([1,5])=[3,1];
nii.original.hdr.dime.dim1([1,5])=[3,1];
movefile(fullfile(NII_dir,sprintf('%s.nii',prefix_nii)),fullfile(NII_dir,sprintf('%s_backup.nii',prefix_nii))); % backup the original
save_nii(nii,fullfile(NII_dir,sprintf('%s.nii',prefix_nii)));

% proessing
for ii=1:1:numel(output_coordinate)
  ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{ii},{prefix_nii;{'backup'}},save_flg);
  niifiles=GetFiles(NII_dir,'*.nii',{prefix_nii;{'backup'}});
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
  ConvertVOIfrom2To1mm(NII_dir,[prefix_nii,'_',output_coordinate{ii}],save_flg);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    delete(fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
    movefile(fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'_highres.voi']),...
             fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
end
delete(fullfile(NII_dir,sprintf('%s.nii',prefix_nii)));
movefile(fullfile(NII_dir,sprintf('%s_backup.nii',prefix_nii)),fullfile(NII_dir,sprintf('%s.nii',prefix_nii))); % restore the backup

CompactVOIs('.',prefix_nii);

cv_hbtools_BVQX_setup(0);

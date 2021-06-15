% script_convert_atlas_NII2VOIs.m
%
% [usage]
% 1. please copy atlas.nii, atlas.txt, and atlas.info onto the directory where this script is located
% 2. please run this script on the MATLAB shell.
%
%
% Created    : "2018-08-27 14:10:30 ban"
% Last Update: "2019-12-16 15:14:43 ban"

cv_hbtools_BVQX_setup(1);

% some constants
NII_dir='.';
labels={};

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii='atlas';
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

% proessing
for ii=1:1:numel(output_coordinate)
  ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{ii},prefix_nii,save_flg);
  niifiles=GetFiles(NII_dir,'*.nii',prefix_nii);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
end
CompactVOIs('.',prefix_nii);

cv_hbtools_BVQX_setup(0);

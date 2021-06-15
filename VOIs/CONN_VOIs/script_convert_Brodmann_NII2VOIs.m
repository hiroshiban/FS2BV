% script_convert_Brodmann_NII2VOIs.m
%
% [usage]
% 1. please copy BA.hdr and BA.img onto the directory where this script is located
% 2. please convert BA.hdr/img files into a nifti file, BA.nii using mricron
% 2. please run this script on the MATLAB shell.
%
%
% Created    : "2018-08-27 11:37:32 ban"
% Last Update: "2018-08-27 16:09:09 ban"

cv_hbtools_BVQX_setup(1);

% some constants
NII_dir='.';
labels={};

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii={'BA';{'atlas','networks'}};
save_flg=1;

% loading ROI label information
fid=fopen('BA.txt','r');
roicounter=0;
tmproi=fgetl(fid);
while ischar(tmproi)
  C=textscan(tmproi,'%d %s');
  roicounter=roicounter+1;
  BA_idx{roicounter}=sprintf('%02d',C{1}(1,:)); %#ok
  BA_hemi{roicounter}=C{2}{1}(1,2:end-2); %#ok
  hemi_idx=strfind(tmproi,C{2}{1}(1,:));
  BA_name{roicounter}=tmproi(hemi_idx+numel(C{2}{1}(1,:))+1:end); %#ok
  tmproi=fgetl(fid);
end
fclose(fid);

[tmp,roi_order]=sort(BA_idx);
BA_idx=BA_idx(roi_order);
BA_hemi=BA_hemi(roi_order);
BA_name=BA_name(roi_order);
nROIs=length(BA_name);

% set ROI colors
lut=[203 110 108; 203 81 94; 203 150 150; 62 62 203; 0 203 0; 162 0 55; 162 55 14;
0 94 92; 203 203 3; 162 0 57; 162 81 203; 203 203 108; 199 199 3; 185 148 39;
203 203 0; 14 176 67; 203 143 203; 0 152 64; 55 55 203; 150 203 203; 141 36 43;
152 69 48; 94 69 55; 39 0 83; 81 0 176; 76 48 127; 0 136 203; 94 18 76;
94 14 14; 94 0 27; 122 203 67; 73 73 201; 203 102 0; 0 203 203; 129 102 0;
203 53 0; 192 24 24; 0 201 203; 108 0 162; 0 201 0; 134 53 0; 203 101 254;
203 0 71; 254 254 5; 169 67 0; 0 250 0; 136 0 203; 254   101   118];

% set labels for ConvertNiftiRoi2BVvoi_Labels
for ii=1:1:nROIs, labels(ii,:)={roi_order(ii),sprintf('BA %s: %s %sH',BA_idx{ii},BA_name{ii},BA_hemi{ii}),lut(str2num(BA_idx{ii}),:)}; end %#ok

% proessing
for ii=1:1:numel(output_coordinate)
  ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{ii},prefix_nii,save_flg);
  niifiles=GetFiles(NII_dir,'*.nii',prefix_nii);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
  ConvertVOIfrom2To1mm(NII_dir,[prefix_nii{1},'_',output_coordinate{ii}],save_flg);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    delete(fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
    movefile(fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'_highres.voi']),...
             fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
end
CompactVOIs('.',prefix_nii);

cv_hbtools_BVQX_setup(0);

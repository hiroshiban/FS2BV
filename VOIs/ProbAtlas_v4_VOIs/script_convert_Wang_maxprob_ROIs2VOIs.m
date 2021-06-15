% script_convert_Wang_maxprob_ROIs2VOIs.m
%
%
% Created    : "2017-08-29 09:23:57 ban"
% Last Update: "2017-09-04 16:20:19 ban"

cv_hbtools_BVQX_setup(1);

% some constants
%NII_dir='./subj_vol_all';
NII_dir='./probability_maps';
labels={};

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii='maxprob_*';
save_flg=1;

labels={1,'V1v';
        2,'V1d';
        3,'V2v';
        4,'V2d';
        5,'V3v';
        6,'V3d';
        7,'hV4';
        8,'VO1';
        9,'VO2';
        10,'PHC1';
        11,'PHC2';
        12,'MST';
        13,'hMT';
        14,'LO2';
        15,'LO1';
        16,'V3b';
        17,'V3a';
        18,'IPS0';
        19,'IPS1';
        20,'IPS2';
        21,'IPS3';
        22,'IPS4';
        23,'IPS5';
        24,'SPL1';
        25,'FEF'};

% proessing
for ii=1:1:numel(output_coordinate)
  ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{ii},prefix_nii,save_flg);
  niifiles=GetFiles(NII_dir,'*.nii.gz',prefix_nii);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
end
ColoringVOIsRainbow(NII_dir);
CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

% script_convert_BrodmannNII2VOIs.m
%
%
% Created    : "2017-09-01 18:49:14 ban"
% Last Update: "2018-08-27 15:55:00 ban"

cv_hbtools_BVQX_setup(1);

% some constants
NII_dir='.';
labels={};

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii='brodmann*';
save_flg=1;

lut=[203	110	108; 203	81	94; 203	150	150; 62	62	203; 0	203	0; 162	0	55; 162	55	14;
0	94	92; 203	203	3; 162	0	57; 162	81	203; 203	203	108; 199	199	3; 185	148	39;
203	203	0; 14	176	67; 203	143	203; 0	152	64; 55	55	203; 150	203	203; 141	36	43;
152	69	48; 94	69	55; 39	0	83; 81	0	176; 76	48	127; 0	136	203; 94	18	76;
94	14	14; 94	0	27; 122	203	67; 73	73	201; 203	102	0; 0	203	203; 129	102	0;
203	53	0; 192	24	24; 0	201	203; 108	0	162; 0	201	0; 134	53	0; 203	101	254;
203	0	71; 254	254	5; 169	67	0; 0	250	0; 136	0	203; 254   101   118];

for ii=1:1:48
  labels(ii,:)={ii,sprintf('Brodmann Area %02d',ii),lut(ii,:)}; %#ok
end

% proessing
for ii=1:1:numel(output_coordinate)
  ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{ii},prefix_nii,save_flg);
  niifiles=GetFiles(NII_dir,'*.nii.gz',prefix_nii);
  for jj=1:1:length(niifiles)
    [niipath,niifname]=fileparts(niifiles{jj});
    movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{ii}),'.voi']));
  end
end
CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

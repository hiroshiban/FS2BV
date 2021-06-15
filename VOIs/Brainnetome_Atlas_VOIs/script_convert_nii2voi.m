% script_convert_nii2voi.m
%
%
% Created    : "2017-08-29 10:05:30 ban"
% Last Update: "2017-09-12 17:52:23 ban"

cv_hbtools_BVQX_setup(1);

% some constants
NII_dir='.';
input_coordinate='MNI';
output_coordinate={'TAL','MNI'};
save_flg=1;

% get ROI lookuptables
colorlut=GetFiles(pwd,'*_LUT.txt');

% processing
for ii=1:1:length(colorlut)
  [lutpath,lutfname]=fileparts(colorlut{ii});
  lutfile=relativepath(colorlut{ii});
  if strcmpi(lutfile(end),filesep()), lutfile=lutfile(1:end-1); end % to omit \ dilimiter at the end of the path
  lut=readFreeSurferColorLUT(lutfile);

  for rr=1:1:length(lut)
    labels(rr,:)={lut{rr}.ID,lut{rr}.name,lut{rr}.RGB};
  end
  prefix_nii=[strrep(lutfname,'_210_LUT',''),'_246'];

  for mm=1:1:length(output_coordinate)
    ConvertNiftiRoi2BVvoi_Labels(NII_dir,labels,input_coordinate,output_coordinate{mm},prefix_nii,save_flg);
    niifiles=GetFiles(NII_dir,'*.nii.gz',prefix_nii);
    for jj=1:1:length(niifiles)
      [niipath,niifname]=fileparts(niifiles{jj});
      movefile(fullfile(niipath,[niifname,'.voi']),fullfile(niipath,[niifname,'_',sprintf('%s',output_coordinate{mm}),'.voi']));
    end
  end
end

CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

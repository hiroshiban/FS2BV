% script_convert_BAAOAnii2VOI.m
%
% usage: function [voi,voi_voxels,nii]=ConvertSPMroi2BVvoi(NII_dir,thres,input_coordinate,output_coordinate,voi_color,prefix_nii,save_flg)
%
%
% Created    : "2019-06-05 10:20:56 ban"
% Last Update: "2019-06-05 14:48:56 ban"

cv_hbtools_BVQX_setup(1);

% constants
input_coordinate='MNI';
output_coordinate={'MNI','TAL'};
prefix_nii='pm';
thres=[0,10];

% processing
for cc=1:1:length(output_coordinate)
  ConvertSPMroi2BVvoi('volume',thres,input_coordinate,output_coordinate{cc},[255,0,0],prefix_nii,1);
  voifiles=GetFiles(fullfile(pwd,'volume'),'*.nii.voi');
  for ii=1:1:length(voifiles)
    [voipath,voifname,voiext]=fileparts(voifiles{ii});
    if ~isempty(strfind(voifname,output_coordinate{cc})), continue; end
    movefile(voifiles{ii},fullfile(voipath,[voifname,'_',output_coordinate{cc},voiext]));
  end
end

% generating 1mm VOIs
ConvertVOIfrom2To1mm('volume','',1);
voifiles=GetFiles(fullfile(pwd,'volume'),'*_highres.voi');
for nn=1:1:length(voifiles)
  movefile(voifiles{nn},strrep(voifiles{nn},'_highres','_1mm'));
end

% renaming
voifiles=GetFiles(fullfile(pwd,'volume'),'*.voi',{{''};{'1mm'}});
for nn=1:1:length(voifiles)
  movefile(voifiles{nn},strrep(voifiles{nn},'.voi','_2mm.voi'));
end

% reduce file sizes
CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

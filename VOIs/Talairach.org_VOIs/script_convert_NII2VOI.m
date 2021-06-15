% script_convert_NII2VOI.m
%
% This simple script converts brain regions defined in Talairach coordinates and saved as NII file
% into BrainVoyager VOI file(s).
%
% [note]
% The original NIFTI-Talairach-coords file (talairach.nii) was obtained from http://www.talairach.org/
% usage 01: ConvertSPMroi2BVvoi(NII_dir,thres,input_coordinate,output_coordinate,voi_color,prefix_nii,save_flg)
% usage 02: CombineVOIs(VOI_dir,save_name,prefix_voi)
% usage 03: ColoringVOIsRainbow(VOI_dir,prefix_voi)
%
%
% Created    : "2016-11-14 09:35:25 ban"
% Last Update: "2018-09-03 12:48:32 ban"

cv_hbtools_BVQX_setup(1);

% constants
NII_dir='.';
nvois=1105;
input_coordinate='TAL';
output_coordinate='TAL';
voi_color='';
prefix_nii='talairach';
save_flg=1;
save_dir='individuals';

if ~exist(save_dir,'dir'), mkdir(save_dir); end

% load nii file and get ROI data
niifile=GetFiles(NII_dir,'*.nii',prefix_nii);
fid=fopen(niifile{1},'r');
dummy=fgetl(fid); % omit header
voi_names=cell(nvois,1);
for ii=1:1:nvois, voi_names{ii}=fgetl(fid); end % read ROI names
fclose(fid);

% processing
for ii=1:1:nvois
  % converting NII ROI to VOI
  fprintf('ROI %04d, ',ii);
  thres=[ii-1,ii];
  ConvertSPMroi2BVvoi(NII_dir,thres,input_coordinate,output_coordinate,voi_color,prefix_nii,save_flg);

  % formating VOI name
  save_voi_name=[strrep(strrep(voi_names{ii},' ','_'),'*',''),'.voi'];
  save_voi_name=strrep(strrep(strrep(save_voi_name,'....voi','.voi'),'...voi','.voi'),'..voi','.voi');
  save_voi_name=strrep(strrep(strrep(save_voi_name,'....',''),'...',''),'..','');
  if strcmp(save_voi_name(1),'.'), save_voi_name=save_voi_name(1,2:end); end

  % saving
  movefile(strrep(niifile{1},'.nii','.voi'),fullfile(save_dir,save_voi_name));

  % change VOI name in the saved VOI file
  voi=BVQXfile(fullfile(save_dir,save_voi_name));
  voi.VOI(1).Name=strrep(save_voi_name,'.voi','');
  voi.Save();
  voi.ClearObject(); clear voi;
end

% combine all the generated VOIs into a file, "talairach.voi"
CombineVOIs(save_dir,'talairach.voi');
movefile(fullfile(save_dir,'talairach.voi'),'.');
ColoringVOIsRainbow('.','talairach');
CompactVOIs('.');
AddPrefix2Filename('.','*.voi','','','_TAL');

cv_hbtools_BVQX_setup(0);

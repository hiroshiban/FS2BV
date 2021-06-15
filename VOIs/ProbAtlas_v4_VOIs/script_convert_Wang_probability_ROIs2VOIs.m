% script_convert_Wang_probability_ROIs2VOIs.m
%
%
% Created    : "2017-08-28 14:20:44 ban"
% Last Update: "2017-09-04 16:20:31 ban"

cv_hbtools_BVQX_setup(1);

% some constants
%NII_dir='./subj_vol_all';
NII_dir='./probability_maps';
thres=[0,25,50,75];

input_coordinate='MNI';
output_coordinate={'TAL','MNI'};

prefix_nii='perc_';
save_flg=1;

ROIs={'V1v','V1d','V2v','V2d','V3v','V3d','hV4','VO1','VO2',...
      'PHC1','PHC2','MST','hMT','LO2','LO1','V3b','V3a','IPS0',...
      'IPS1','IPS2','IPS3','IPS4','IPS5','SPL1','FEF'};

% proessing
for ii=1:1:numel(thres)
  for jj=1:1:length(output_coordinate)
    niifiles=GetFiles(NII_dir,'perc_*.nii.gz');
    for vv=1:1:length(niifiles)
      [niipath,niifname]=fileparts(niifiles{vv});
      sidx=strfind(niifname,'_roi');
      sidx=sidx+4;
      eidx=strfind(niifname,'.nii');
      eidx=eidx-4;
      if ~isempty(strfind(niifname,'_lh'))
        hemi='_lh';
      else
        hemi='_rh';
      end
      movefile(niifiles{vv},fullfile(niipath,[ROIs{str2num(niifname(sidx:eidx))},hemi,'.nii.gz']));
      gunzip(fullfile(niipath,[ROIs{str2num(niifname(sidx:eidx))},hemi,'.nii.gz']));
      delete(fullfile(niipath,[ROIs{str2num(niifname(sidx:eidx))},hemi,'.nii.gz']));
    end
    ConvertSPMroi2BVvoi(NII_dir,[thres(ii),100],input_coordinate,output_coordinate{jj},[],{'_lh','_rh'},save_flg);
    CombineVOIs(NII_dir,sprintf('WangProbAtlas_thr%02d_%s.voi',thres(ii),output_coordinate{jj}),{'_lh','_rh'});
    delete(fullfile(NII_dir,'*_lh*.voi'));
    delete(fullfile(NII_dir,'*_rh*.voi'));
  end
end
ColoringVOIsRainbow(NII_dir);
CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

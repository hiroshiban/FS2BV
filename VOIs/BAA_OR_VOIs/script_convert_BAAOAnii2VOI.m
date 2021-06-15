% script_convert_BAAOAnii2VOI.m
%
% usage: function [voi,voi_voxels,nii]=ConvertSPMroi2BVvoi(NII_dir,thres,input_coordinate,output_coordinate,voi_color,prefix_nii,save_flg)
%
%
% Created    : "2019-06-05 10:20:56 ban"
% Last Update: "2019-12-16 15:14:36 ban"

cv_hbtools_BVQX_setup(1);

% constants (note: all ROI names are from documentations of the atlas)
NII_dir={'face','object','scene'};
ROIs{1}={'rOFA','lOFA','rpFus','lpFus','raFus','laFus','rpcSTS','lpcSTS','rpSTS','lpSTS','raSTS','laSTS'};
ROIs{2}={'rLO','lLO','rpFs','lpFs'};
ROIs{3}={'rTOS','lTOS','rRSC','lRSC','rPPA','lPPA'};

input_coordinate={'MNI'};
output_coordinate={'MNI','TAL'};
prefix_nii='-thr10-2mm';

% processing
cdir=fileparts(mfilename('fullpath'));
for ii=1:1:length(NII_dir)
  cd(fullfile(NII_dir{ii},'volume'));

  niifile=GetFiles(pwd,[prefix_nii,'*.nii.gz']);
  [niipath,niifname]=fileparts(niifile{1});

  voi_color=int8(round(255.*generate_colormap_rainbow(length(ROIs{ii}))));
  for cc=1:1:length(output_coordinate)
    for rr=1:1:length(ROIs{ii})
      thres=[rr-1+0.001,rr];
      ConvertSPMroi2BVvoi('.',thres,input_coordinate{1},output_coordinate{cc},voi_color(rr,:),prefix_nii,1);
      movefile(strrep(niifile{1},'.nii.gz','.nii.voi'),sprintf('%s_%s.voi',ROIs{ii}{rr},output_coordinate{cc}));
    end
    for rr=1:1:length(ROIs{ii})
      voi=BVQXfile(sprintf('%s_%s.voi',ROIs{ii}{rr},output_coordinate{cc}));
      voi.VOI(1).Name=ROIs{ii}{rr};
      voi.Save();
    end
    CombineVOIs('.',[niifname,sprintf('_%s.voi',output_coordinate{cc})],{{''};{prefix_nii}});
    for rr=1:1:length(ROIs{ii}), delete(fullfile(pwd,sprintf('%s_%s.voi',ROIs{ii}{rr},output_coordinate{cc}))); end
  end
  cd(cdir);
end % for ii=1:1:length(NII_dir)

% generating 1mm VOIs
cdir=fileparts(mfilename('fullpath'));
for ii=1:1:length(NII_dir)
  cd(fullfile(NII_dir{ii},'volume'));
  ConvertVOIfrom2To1mm('.','',1);
  voifiles=GetFiles(pwd,['_highres.voi']);
  for nn=1:1:length(voifiles), movefile(voifiles{nn},strrep(strrep(voifiles{nn},'_highres',''),'2mm','1mm')); end
  cd(cdir);
end

% reduce file sizes
CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

% script_script_coloring_HCP_MMP1_VOIs.m
%
% This script converts FSL HCP MMP1 NIFTI atlas (ROIs, e.g. HCP-MMP_1mm.nii.gz) to BrainVoyager VOIs.
% Then, the script paints the generated HCP MMP1 VOIs with the colors defined in {lh|rh}.HCP-MMP1.annot
% files diestributed in
% https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446
%
% see also:
% https://osf.io/azup8
%
%
% Created    : "2024-03-15 17:46:39 ban"
% Last Update: "2024-07-09 23:36:47 ban"

% loading BVQX_hbtools
cv_hbtools_BVQX_setup(1);

% converting FLS ROI (*.nii) to BrainVoyager VOIs.
ConvertFSLroi2BVvoi('.','MNI','MNI','*cortices',1);
ConvertFSLroi2BVvoi('.','MNI','MNI',{'';'*cortices'},1);

% VOI coloring
% from here the script paints the HCP MMP1 VOIs generated from nifti files (e.g. HCP-MMP_1mm.nii.gz)

% some constants
%VOIs={'HCP-MMP1/HCP-MMP1_1mm.nii.voi','HCP-MMP1/HCP-MMP1_2mm.nii.voi',...
%      'HCP-MMP1_cortices/HCP-MMP1_cortices_1mm.nii.voi','HCP-MMP1_cortices/HCP-MMP1_cortices_2mm.nii.voi'};
VOIs=wildcardsearch('.','*.voi');

% loading colortable_lh and colortable_rh
load(fullfile(fileparts(mfilename('fullpath')),'HCP_MMP1_annotation','FreeSurfer_HCP_MMP_ColorLUT.mat'));

% processing
for ii=1:1:length(VOIs)
  fprintf('processing: %s...',VOIs{ii});
  voi=BVQXfile(fullfile(fileparts(mfilename('fullpath')),VOIs{ii}));
  for vv=1:1:voi.NrOfVOIs
    for rr=1:1:colortable_lh.numEntries
      if strcmp(voi.VOI(vv).Name,strcat(strrep(colortable_lh.struct_names{rr}(3:end),'_ROI',''),'_L'))
        voi.VOI(vv).Color=colortable_lh.table(rr,1:3);
      end
    end

    for rr=1:1:colortable_rh.numEntries
      if strcmp(voi.VOI(vv).Name,strcat(strrep(colortable_rh.struct_names{rr}(3:end),'_ROI',''),'_R'))
        voi.VOI(vv).Color=colortable_rh.table(rr,1:3);
      end
    end
  end
  voi.SaveAs(fullfile(fileparts(mfilename('fullpath')),VOIs{ii}));
  fprintf('done.\n');
end

% unloading BVQX_hbtools
cv_hbtools_BVQX_setup(0);

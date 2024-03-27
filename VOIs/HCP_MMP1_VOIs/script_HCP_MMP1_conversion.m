% script_HCP_MMP1_conversion.m
%
% a simple script to read HCP MMP1 ROI XMLs and the corresponding *.nii.gz files,
% and convert them to BrainVoyager VOI files.
%
% [ref]
% Glasser, M. F., Coalson, T. S., Robinson, E. C., Hacker, C. D., Harwell, J., Yacoub, E., et al. (2016).
% A multi-modal parcellation of human cerebral cortex. Nature, 1-11.
% http://doi.org/10.1038/nature18933
%
%
% Created    : "2024-03-27 21:12:31 ban"
% Last Update: "2024-03-27 21:27:47 ban"

cv_hbtools_BVQX_setup(1);

ConvertFSLroi2BVvoi('.','MNI','MNI','*cortices',1);
ConvertFSLroi2BVvoi('.','MNI','MNI',{'';'*cortices'},1);

CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

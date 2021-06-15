% script_convert_AALnii2VOI.m
%
%
% Created    : "2017-08-28 14:08:57 ban"
% Last Update: "2017-08-31 17:48:01 ban"

% IMPORTANT NOTE
% To reslice the AAL.nii from 2mm to 1mm cubic space, please run
% >> reslice_nii('AAL.nii','AAL_1mm.nii',1,'','',2);
% reslice_nii is in BVQX_hbtools/nifti_tools.

cv_hbtools_BVQX_setup(1);

ConvertAALroi2BVvoi('.','MNI','TAL','',1);
movefile('AAL2.voi','AAL2_TAL.voi');
movefile('AAL2_1mm.voi','AAL2_TAL_1mm.voi');

ConvertAALroi2BVvoi('.','MNI','MNI','',1);
movefile('AAL2.voi','AAL2_MNI.voi');
movefile('AAL2_1mm.voi','AAL2_MNI_1mm.voi');

CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

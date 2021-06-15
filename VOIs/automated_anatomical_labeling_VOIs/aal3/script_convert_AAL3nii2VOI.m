% script_convert_AALnii2VOI.m
%
%
% Created    : "2021-06-09 10:38:41 ban"
% Last Update: "2021-06-09 10:39:44 ban"

cv_hbtools_BVQX_setup(1);

ConvertAALroi2BVvoi('.','MNI','TAL','',1);
movefile('AAL3v1.voi','AAL3v1_TAL.voi');
movefile('AAL3v1_1mm.voi','AAL3v1_TAL_1mm.voi');
movefile('ROI_MNI_V7.voi','ROI_MNI_V7_TAL.voi');
movefile('ROI_MNI_V7_1mm.voi','ROI_MNI_V7_TAL_1mm.voi');

ConvertAALroi2BVvoi('.','MNI','MNI','',1);
movefile('AAL3v1.voi','AAL3v1_MNI.voi');
movefile('AAL3v1_1mm.voi','AAL3v1_MNI_1mm.voi');
movefile('ROI_MNI_V7.voi','ROI_MNI_V7_MNI.voi');
movefile('ROI_MNI_V7_1mm.voi','ROI_MNI_V7_MNI_1mm.voi');

CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

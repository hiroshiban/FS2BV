% script_convert_TAL2MNI_VOI.m
%
%
% Created    : "2017-08-27 14:50:28 ban"
% Last Update: "2018-09-03 12:48:27 ban"

cv_hbtools_BVQX_setup(1);

voi=BVQXfile('shen_parcels_TAL.voi');
for ii=1:1:voi.NrOfVOIs, voi.VOI(ii).Voxels=tal2mni(voi.VOI(ii).Voxels); end
voi.ReferenceSpace='MNI';
voi.SaveAs('shen_parcels_MNI.voi');
voi.ClearObject(); clear voi;

cv_hbtools_BVQX_setup(0);

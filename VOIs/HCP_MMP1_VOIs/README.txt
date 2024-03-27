****************************************
README.txt on atlases

Created    : "2024-03-27 21:14:11 ban"
Last Update: "2024-03-27 21:17:18 ban"
****************************************

The BrainVoyager VOI files in this directory are ones extracted from HCP-MMP1 atlases
converted from *.nii.gzip to *.voi (the MNI coordinate), using ConvertFSLroi2BVvoi.m
function in BVQX_hbtools.

Please note that some VOI files contains more than 255 VOIs in one *.voi file, while the
current version of BrainVoyager can support up to 255 VOIs. Therefore, you can not use
some VOI files as they are but need to extract/delete what you need before using them in
your analyses.

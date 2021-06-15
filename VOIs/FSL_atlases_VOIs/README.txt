****************************************
README.txt on atlases

Created    : "2015-12-25 16:22:36 ban"
Last Update: "2015-12-26 15:45:58 ban"
****************************************

The BrainVoyager VOI files in this directory are ones extracted from FSL atlases whose
coordinates are converted from MNI to Talairach, using ConvertFSLroi2BVvoi.m function in
BVQX_hbtools.

Please note that some VOI files contains more than 255 VOIs in one *.voi file, while the
current version of BrainVoyager can support up to 255 VOIs. Therefore, you can not use
some VOI files as they are but need to extract/delete some before using them in your analyses.

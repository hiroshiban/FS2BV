% script_FSL_ROI_conversion.m
%
% a simple script to read FSL ROI XMLs and the corresponding *.nii.gz files,
% and convert them to BrainVoyager VOI files.
%
%
% Created    : "2015-12-26 15:44:06 ban"
% Last Update: "2017-09-05 12:57:19 ban"

cv_hbtools_BVQX_setup(1);

% some constants
XML_dir={'./TAL','./MNI'};
input_coordinate={'MNI','MNI'};
output_coordinate={'TAL','MNI'};
prefix_xml='';
save_flg=1;

% processing
for ii=1:1:length(XML_dir)
  ConvertFSLroi2BVvoi(XML_dir{ii},input_coordinate{ii},output_coordinate{ii},prefix_xml,save_flg);
  AddPrefix2Filename(XML_dir{ii},'*.voi','','',sprintf('_%s',strrep(XML_dir,'./','')));
end

CompactVOIs('.');

cv_hbtools_BVQX_setup(0);

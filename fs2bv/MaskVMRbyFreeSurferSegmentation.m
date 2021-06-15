function MaskVMRbyFreeSurferSegmentation(vmr_file,freesurfer_wm_seg_file,flip_flg,vmr_framingcube)

% Masks the white matter region of the input BrainVoayger's native space VMR by FreeSurfer's segmentation (wm.seg.mgz)
% function MaskVMRbyFreeSurferSegmentation(vmr_file,:freesurfer_wm_seg_file,:flip_flg,:vmr_framingcube)
% (: is optional)
%
% This function reads native(or ACPC)-space VMR file and masks its white-matter regions by
% the FreeSurfer's segmentation result (wm.seg.mgz in the $FREESURFER_SUBJECT_DIR/(subj)/mri directory).
% The masked VMR can be used for fine-scale gray/whiete-matter segmentation and cortical reconstruction.
% Here, note that a masked VMR generally provides better segmentation result, but it depends on
% the image quality (in some cases, BrainVoyager segmentation alone gives more accurate results).
% Please be careful in use.
%
% [about the difference between MaskVMRbyFreeSurferSegmentation and MaskVMRbyFreeSurferSegmentation_ribbon]
% MaskVMRbyFreeSurferSegmentation       : for general masking purposes. Any *.mgz segmentation result can be used as a
%                                         mask (by default, the parameters are tuned to process wm.seg.mgz as a mask)
% MaskVMRbyFreeSurferSegmentation_ribbon: Specific for applying a mask using the white (and gray) matter segmentation
%                                         result in ribbon.mgz.
%                                         Generally, for surface reconstructions, MaskVMRbyFreeSurferSegmentation_ribbon
%                                         gives the better results.
%
% [NOTE on importing ROIs, segmentation, surfaces etc defined by the other software, such as FSL or FreeSurfer into BrainVoyager]
% To import the ROIs defined outside BrainVoyager, I prepared the functions below.
% ConvertNiftiRoi2BVvoi_ProbThres : Converts NII-format ROI probability map to BrainVoayer VOIs
%                                   with thresholding the map values
% ConvertNiftiRoi2BVvoi_Labels    : Converts NII-format ROI probability map to BrainVoayer VOIs
%                                   using the label lookuptable corresponding to the map ID
% ExtractFSLroi            : Extracts specific value(s) from NII based on XML database
% ExtractFSLroiDirect      : Extracts specific value(s) from NII directly for ROI generations
% ConvertFSLroi2BVvoi      : Converts FSL NII ROIs to BrainVoyager VOIs
% ConvertSPMroi2BVvoi      : Converts SPM NII ROIs to BrainVoyager VOIs
% ConvertsAALroi2BVvoi     : Converts SPM AAL antomical tempolate (NII) ROIs to BrainVoyager VOIs
% ConvertFreeSurferAnnotation2BVpoi   : Converts FreeSurfer surface annotations to BrainVoyager POIs
% ConvertFreeSurferParcellation2BVvoi : Converts FreeSurfer MGZ parcellations to BrainVoyager VOIs
% ConvertFreeSurferMGZ2VMR : Converts FreeSurer MGZ T1/ROI files to BrainVoayer VMRs
% ConvertFreeSurferRibbon2BL2VMR : Converts FreeSurfer ribbon.mgz to BrainVoyager *_{LH|RH}_BL2.vmr
% ConvertFreeSurferSurf2SRF: Converts FreeSurer surface files to BrainVoayer SRFs
% ImportFreeSurfer2BrainVoyager : Imports FreeSurfer-processed files into BrainVoyager
% ConvertBVpoi2FreeSurferAnnotation   : Converts BrainVoyager POIs to FreeSurfer Annotation files. We can further
%                                       generate label or volume ROI files from the generated annotation files using
%                                       FreeSurfer commands.
%
% Some of the VOIs defined in TAL/MNI spaces are already stored in ~/BVQX_hbtools/VOIs.
% To find the VOIs in which specific XYZ coordinates in TAL/MNI space, please use the function below,
% GetAreaNameFromAtlasVOI               : Returns area candidates, in which the input XYZ coordinate(s)
%                                        is(are) belonging to, based on the pre-defined VOI atlases.
%
% [input]
% vmr_file : BrainVoyager VMR file (in the native or ACPC space, not TAL or MNI)
%            The file should be specified with a relative path format in which
%            the location where this function is called is the origin of the path.
%            !!!IMPORTANT NOTE!!!
%            Currently, only 1mm cubic VMR can be accepted. Please be careful.
% freesurfer_wm_seg_file : (optional) FreeSurfer's segmentation result file,
%            that is originally stored as $FREESURFER_SUBJECT_DIR/(subj)/mri/wm.seg.mgz
%            after running the recon-all function. After copying the file somewhere,
%            please specify its location by this variable with a relative path format.
%            e.g. freesurfer_wm_seg_file='../3d/wm.seg.mgz';
%            if not specified, freesurfer_wm_seg_file='./wm.seg.mgz' is used by default.
% flip_flg : (optional) how to flip the FreeSurfer coordinate so that it matches with the
%            BrainVoyager VMR coordinate.
%            if 0, the FreeSurfer-processed anatomy is flipped along x- and z- axes.
%            if non-zero, the FreeSurfer-processed anatomy is flipped only along x-axis.
%            1 by default. Generally, if we need to set this value 0, the VMR would be
%            flipped. Please be careful.
% vmr_framingcube : (optional) Framing cube dimension (pixels) of the VMR to be generated.
%            Generally this value is one of 256, 386, and 512. 256 by default (1mm VMR's defalt in BV).
%
% [output]
% no output variable
% The input VMR file is overwritten by the masked VMR file and the original
% file is backuped with '_backup' prefix in the same directory.
%
% [dependency]
% UNIX command line tools such as gzip are required to use this function.
% Especially, if you want to use this function on Windows, please install UNIX command emulator,
% Cygwin or MSYS2, and set an envitonmental path to the tools in advance.
%
%
% Created    : "2017-07-08 16:56:26 ban"
% Last Update: "2018-07-25 09:36:18 ban"

% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(freesurfer_wm_seg_file), freesurfer_wm_seg_file='./wm.seg.mgz'; end
if nargin<3 || isempty(flip_flg), flip_flg=1; end
if nargin<4 || isempty(vmr_framingcube), vmr_framingcube=256; end

if ~exist(fullfile(pwd,vmr_file),'file')
  error('The input VMR file not found. check the input variable.');
end

if ~exist(fullfile(pwd,freesurfer_wm_seg_file),'file')
  error('The input FreeSurfer segmentation file not found. check the input variable.');
end

% check whether the input vmr xyz dimensions are isometric
vmr=BVQXfile(fullfile(pwd,vmr_file));
if vmr.DimX~=vmr.DimY || vmr.DimY~=vmr.DimZ
  error('VMR XYZ voxel dimensions mismatched. only a cubic VMR is accepted. check the data');
end

% set a path to the tools for handling FreeSurfer files on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

[vmrpath,vmrfname,vmrext]=fileparts(fullfile(pwd,vmr_file));
[segpath,segfname,segext]=fileparts(fullfile(pwd,freesurfer_wm_seg_file));

fprintf('Target BrainVoyager VMR            : %s%s\n',vmrfname,vmrext);
fprintf('Target FreeSurfer segmentation file: %s%s\n',segfname,segext);
fprintf('\n');

%% processing FreeSurfer wm.seg.mgz file

%% === 1. converting wm.seg.mgz to wm.seg.vmr ===
fprintf('Converting FreeSurfer segmentation file to BrainVoyager VMR...');
segfile=MRIread(fullfile(pwd,freesurfer_wm_seg_file));
if ~isempty(find(uint8(segfile.volres(1))~=1,1))
  error('currently only 1mm cubic segmentation file is accepted. check the input FreeSurfer segmentation file.');
end

segvmr=BVQXfile('new:vmr');
segvmr.FramingCube=vmr_framingcube;
segvmr.DimX=segfile.height;
segvmr.DimY=segfile.width;
segvmr.DimZ=segfile.depth;
segvmr.NRows=segvmr.DimX;
segvmr.FOVRows=segvmr.DimX;
segvmr.FOVCols=segvmr.DimY;

% note: here, to convert the coordinate from FreeSurfer to BrainVoyager, the volume should be
% 1. adjusted contrast and brightness of the segmented white matter region around 180.
% 2. transformed from the default FreeSurfer coordinate into the BrainVoyager SAG coordinate.

% brightness adjustment
sz=size(segfile.vol);
vimg=segfile.vol(:);
imgidx=find(vimg>50);
subimg=vimg(imgidx); % process only regions with value > 50
subimg=subimg-min(subimg(:)); % adjusting subimg within 0.0-1.0
subimg=subimg./max(subimg(:));
subimg=10.*subimg-5+180; % sdjusting subimg within 175-185
vimg(imgidx)=subimg;
vimg=uint8(reshape(vimg,sz));

% transformation
try
  if ~flip_flg
    segvmr.VMRData=flip(flip(permute(vimg,[3,1,2]),1),3);
  else
    segvmr.VMRData=flip(permute(vimg,[3,1,2]),1);
  end
catch
  if ~flip_flg
    segvmr.VMRData=flipdim(flipdim(permute(vimg,[3,1,2]),1),3);
  else
    segvmr.VMRData=flipdim(permute(vimg,[3,1,2]),1);
  end
end

% saving the result
segvmr.SaveAs(fullfile(segpath,[segfname,'.vmr']));
fprintf('done.\n');

%% === 2. coregistration of FreeSurfer segmentation VMR to the input BrainVoyager VMR.
fprintf('\n');
fprintf('Please follow the instructions below.\n\n');
fprintf('1. Coregist FreeSurfer segmentation %s.vmr file\n',segfname);
fprintf('   to the input VMR using BrainVoyager\n');
fprintf('   and save the transformation matrix file as.\n');
fprintf('   %s-TO-%s.trf.\n',segfname,vmrfname);
%fprintf('2. Using the transformation matrix file above,\n');
%fprintf('   convert the FreeSurfer segmentation %s.vmr file\n',segfname);
%fprintf('   and save it as %s_TRF.vmr.\n',segfname);
fprintf('Then, press F5 (or type dbcont) to proceed\n');
keyboard;

fprintf('applying spatial transformation from FreeSurfer native MGZ to BrainVoyager VMR space...');
trf=BVQXfile(fullfile(segpath,sprintf('%s-TO-%s.trf',segfname,vmrfname)));
tmpvmr=segvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
tmpvmr.SaveAs(fullfile(segpath,[segfname,'_TRF.vmr']));
trf.ClearObject(); clear trf;
tmpvmr.ClearObject(); clear tmpvmr;
fprintf('done.\n');

%% === 3. mask the white matter region of the input VMR by the FreeSurfer segmentation file
fprintf('masking %s%s...',vmrfname,vmrext);
vmr=BVQXfile(fullfile(pwd,vmr_file));
mask=BVQXfile(fullfile(segpath,[segfname,'_TRF.vmr']));

% the second brightness adjustment. since the image intensities at the edges
% are interpolated after the transformation above, we have to adjust image
% contrast and brightness again here.
sz=size(mask.VMRData);
vimg=mask.VMRData(:);
imgidx=find(vimg>50);
subimg=vimg(imgidx); % process only regions with value > 50
subimg=subimg-min(subimg(:)); % adjusting subimg within 0.0-1.0
subimg=subimg./max(subimg(:));
subimg=10.*subimg-5+180; % sdjusting subimg within 175-185
vimg(imgidx)=subimg;
vimg=uint8(reshape(vimg,sz));

vmr.SaveAs(fullfile(vmrpath,[vmrfname,'_backup.vmr']));
vmr.VMRData(vimg>170)=vimg(vimg>170);
vmr.SaveAs(fullfile(vmrpath,[vmrfname,vmrext]));
fprintf('done.\n');

% clean up
vmr.ClearObject(); clear vmr;
mask.ClearObject(); clear mask;

% remove a path to the tools for handling FreeSurfer files on MATLAB
rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

return

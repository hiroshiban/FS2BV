function vmr=ConvertFreeSurferRibbon2BL2VMR(ribbon_file,VMR_coregistration_file,ACPC_file,TAL_file,MNI_file,flip_flg,vmr_framingcube)

% Reads a FreeSurfer ribbon.mgz file and converts it to BrainVoyager *{LH|RH}_BL2_*.vmr (gray/white-matter segmentation) files.
% function vmr=ConvertFreeSurferRibbon2BL2VMR(ribbon_file,:VMR_coregistration_file,:ACPC_file,:TAL_file,:MNI_file,:flip_flg,:vmr_framingcube)
% (: is optional)
%
% This function reads a FreeSurfer ribbon.mgz file and converts it to a BrainVoyager VMR file.
% If you just want to convert FreeSurfer MGZ to BrainVoyager VMR format, please set mgz_file alone and run this
% function without any other input variables. If you want to generate the VMR files in the BrainVoyager-native
% subject anatomy space, please additionally set the 'VMR_coregistration_file' variable. If you want to create
% VMRs in ACPC, TAL, or MNI space, please further set 'ACPC_file' or/and 'TAL/MNI_file' respectively. To create
% TAL VMRs, you don't need to set MNI_file, while you don't need to set TAL_file (so please set TAL_file='';)
% to create MNI VMRs.
%
% To create the corresponding T1 anatomy file (T1.mgz.vmr), for instatnce, plesae run the commands below
% >> mgz_file=fullfile(freesurfer_subj_dir,'HB','mri','T1.mgz');
% >> ConvertFreeSurferMGZ2VMR(mgz_file);
%
% For FreeSurfer details, please see
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki
% https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
% https://surfer.nmr.mgh.harvard.edu/fswiki/SubcorticalSegmentation
% https://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011
%
% !!! IMPORTANT NOTE 1 !!!
% It may take relatively long computation time to convert a FreeSurfer MGZ into VMR in TAL space
% using MATLAB and this function. Therefore it would be recommended to create a VMR file in ACPC
% space and convert it into TAL space using the transformation option implemented in the BrainVoyager's
% Volume GUI window.
%
% !!! IMPORTANT NOTE 2 !!!
% After the Talairach transformation, the original volume values, which are assined as regional/tissue
% IDs are not preserved precisely due to spatial interpolations (native and ACPC conversions are fine
% as I have used a nearest neighbor interpolation method). So if you want to preserve the ID accurately,
% please convert the MGZ to ACPC VMR using this function and then convert the ACPC VMR to TAL using
% BrainVoyager GUI window, which allows us to convert VMR to the Talairach space using a nearest neighbor
% sampling method.
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
% MaskVMRbyFreeSurferSegmentation       : for general masking purposes. Any *.mgz segmentation result can be used as a
%                                         mask (by default, the parameters are tuned to process wm.seg.mgz as a mask)
% MaskVMRbyFreeSurferSegmentation_ribbon: Specific for applying a mask using the white (and gray) matter segmentation
%                                         result in ribbon.mgz.
%                                         Generally, for surface reconstructions, MaskVMRbyFreeSurferSegmentation_ribbon
%                                         gives the better results.
% ConvertBVpoi2FreeSurferAnnotation   : Converts BrainVoyager POIs to FreeSurfer Annotation files. We can further
%                                       generate label or volume ROI files from the generated annotation files using
%                                       FreeSurfer commands.
%
% Some of the VOIs defined in TAL/MNI spaces are already stored in ~/BVQX_hbtools/VOIs.
% To find the VOIs in specific XYZs of TAL/MNI coords, please use the function below,
% GetAreaNameFromAtlasVOI               : Returns area candidates, in which the input XYZ coordinate(s)
%                                        is(are) belonging to, based on the pre-defined VOI atlases.
%
% Some *.nii/*.nii.gz/*.img files are provided as ROIs 2mm cubic volume space.
% If you want to convert them to 1mm cubic space, that is compatible with BrainVoayger default anatomical
% volume space, please run the command below
% >> reslice_nii('AAL.nii','AAL_1mm.nii',1,'','',2);
% before using the functions above. 'reslice_nii' is in BVQX_hbtools/nifti_tools.
%
% [example]
% >> ribbon_file='./freesurfer/subjects/DC/mri/ribbon.mgz';
% >> VMR_coregistration_file='./DC/3d/wm.seg-TO-DC_10_dc17_006.3d_final.trf';
% >> ACPC_file='./DC/3d/DC_10_dc17_006.3d_final_ACPC.trf';
% >> vmr=ConvertFreeSurferRibbon2BL2VMR(ribbon_file,VMR_coregistration_file,ACPC_file);
%
% [input]
% mgz_file    : A FreeSurfer MGZ file to be converted to VMR specified with a relative path format,
%               in which the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC/aparc.a2009s+aseg.mgz';
% VMR_coregistration_file : (optional) a transformation matrix file that is generated after coregistering
%               (fs_subj_dir)/mri/wm.seg.mgz(vmr) to the corresponding subject's BrainVoyager VMR
%               file in the native space (generally *_final.vmr in BVQX_hbtools).
%               for the details of coregistration, please see MaskVMRbyFreeSurferSegmentation.m
%               The file should be specified as a relative path format.
% ACPC_file   : (optional) BrainVoyager ACPC transformation file (*_ACPC.trf), specified
%               as a relative path format. If not specified, VMRs are generated in the native
%               space. empty by default.
% TAL_file    : (optional) BrainVoyager TAL transformation file (*_TAL.tal), specified
%               as a relative path format. If not specified, VMRs are generated in the native
%               or ACPC space. To generate a VOI file in TAL space, please also specify ACPC_file.
%               empty by default.
% MNI_file    : (optional) BrainVoyager MNI transformation file (*_TO_MNI_a12.trf), specified
%               as a relative path format. If not specified, VMRs are generated in the native,
%               ACPC, or TAL space. To generate a VMR file in MNI space, please also specify
%               ACPC_file, but no need to specify TAL_file (TAL_file=''; is fine).
%               empty by default.
% flip_flg    : (optional) how to flip the FreeSurfer coordinate so that it matches with the
%               BrainVoyager VMR coordinate.
%               if 0, the FreeSurfer-processed anatomy is flipped along x- and z- axes.
%               if 1, the FreeSurfer-processed anatomy is flipped only along x-axis.
%               1 by default. Generally, if we need to set this value 0, the VMR would be
%               flipped. Please be careful.
% vmr_framingcube : (optional) Framing cube dimension (pixels) of the VMR to be generated.
%               Generally this value is one of 256, 386, and 512. 256 by default (1mm VMR's defalt in BV).
%
% [output]
% vmr         : A BrainVoyager VMR object.
% The generated VMR file is saved in the same directory with the input MGZ file with an extention, '_{LH|RH}_BL2.vmr'.
%
% [dependency]
% 1. freesurfer_matlab_tools
% Matlab tools to read/write FreeSurfer files.
% The tools are copied under ~/freesurfer/matlab/ when you install FreeSurfer.
%
% 2. UNIX command line tools such as gzip are required to use this function.
% Especially, if you want to use this function on Windows, please install UNIX command emulator,
% Cygwin or MSYS2, and set an envitonmental path to the tools in advance.
%
%
% Created    : "2017-09-14 11:08:28 ban"
% Last Update: "2021-06-16 10:02:40 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(VMR_coregistration_file), VMR_coregistration_file=''; end
if nargin<3 || isempty(ACPC_file), ACPC_file=''; end
if nargin<4 || isempty(TAL_file), TAL_file=''; end
if nargin<5 || isempty(MNI_file), MNI_file=''; end
if nargin<6 || isempty(flip_flg), flip_flg=1; end
if nargin<7 || isempty(vmr_framingcube), vmr_framingcube=256; end

% check the input files
if ~exist(fullfile(pwd,ribbon_file),'file'), error('ribbon_file not found. check the input variable.'); end

if isempty(ACPC_file) && ~isempty(TAL_file)
  error('if you want to generate VOI in TAL space, you have to set ACPC_file too.');
end

if ~isempty(VMR_coregistration_file) && ~exist(fullfile(pwd,VMR_coregistration_file),'file')
  error('VMR_coregistration_file not found. check the input variable.');
end

if ~isempty(ACPC_file) && ~exist(fullfile(pwd,ACPC_file),'file')
  error('ACPC_file not found. check the input variable.');
end

if ~isempty(TAL_file) && ~exist(fullfile(pwd,TAL_file),'file')
  error('TAL_file not found. check the input variable.');
end

if ~isempty(MNI_file) && ~exist(fullfile(pwd,MNI_file),'file')
  error('MNI_file not found. check the input variable.');
end

output_coordinate='native';
if ~isempty(VMR_coregistration_file), output_coordinate='native'; end
if ~isempty(ACPC_file), output_coordinate='ACPC'; end
if ~isempty(TAL_file), output_coordinate='TAL'; end
if ~isempty(MNI_file), output_coordinate='MNI'; end

% set a path to the tools for handling FreeSurfer files on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

% display message
fprintf('Target FreeSurfer Ribbon MGZ : %s\n',ribbon_file);
fprintf('Output coordinate            : %s\n',output_coordinate);
fprintf('\n');

%% processing

% loading the input MGZ file
[mgzpath,mgzfname,mgzext]=fileparts(fullfile(pwd,ribbon_file));
segfile=MRIread(fullfile(pwd,ribbon_file));

% initializing and assigining the segmentation values to the VMRs
% generating an output vmr file
hemis={'LH','RH'};
vmr=cell(2,1);
for hh=1:1:length(hemis)
  mgzvmr=BVQXfile('new:vmr');
  mgzvmr.CoordinateSystem=2;
  mgzvmr.ReferenceSpace=1; % native coordinate
  mgzvmr.VoxResVerified=0;
  mgzvmr.FramingCube=vmr_framingcube;
  mgzvmr.DimX=segfile.height;
  mgzvmr.DimY=segfile.width;
  mgzvmr.DimZ=segfile.depth;
  mgzvmr.NRows=mgzvmr.DimX;
  mgzvmr.FOVRows=mgzvmr.DimX;
  mgzvmr.FOVCols=mgzvmr.DimY;

  % assign white and gray matter segmentations as VMRdata
  if hh==1 % left hemisphere
    % extracting the regions from the left hemisphere, here ID=2 is white matter and 3 is gray matter,
    % which are defined in FreeSurferColorLUT.txt
    vol=uint8(zeros(size(segfile.vol))); % for left hemisphere
    vol(segfile.vol==2)=240;
    %vol(segfile.vol==3)=120;%235;
  else % right hemisphere
    % extracting the regions from the left hemisphere, here ID=41 is white matter and 42 is gray matter,
    % which are defined in FreeSurferColorLUT.txt
    vol=uint8(zeros(size(segfile.vol))); % for left hemisphere
    vol(segfile.vol==41)=240;
    %vol(segfile.vol==42)=120;%235;
  end
  mgzvmr.VMRData=uint8(vol);
  clear vol;

  % transformation
  % note: here, to convert the coordinate from FreeSurfer to BrainVoyager, the volume should be
  % 1. adjusted contrast and brightness of the segmented white matter region around 180.
  % 2. transformed from the default FreeSurfer coordinate into the BrainVoyager SAG coordinate.
  try
    if ~flip_flg
      mgzvmr.VMRData=flip(flip(permute(mgzvmr.VMRData,[3,1,2]),1),3);
    else
      mgzvmr.VMRData=flip(permute(mgzvmr.VMRData,[3,1,2]),1);
    end
  catch
    if ~flip_flg
      mgzvmr.VMRData=flipdim(flipdim(permute(mgzvmr.VMRData,[3,1,2]),1),3);
    else
      mgzvmr.VMRData=flipdim(permute(mgzvmr.VMRData,[3,1,2]),1); %#ok
    end
  end

  % transform the space
  if ~strcmpi(output_coordinate,'native') || ~isempty(VMR_coregistration_file)

    % transformation using the result of VMR (wm.seg.vmr) to VMR (*_final.vmr) coregistration
    if ~isempty(VMR_coregistration_file)
    %if strcmpi(output_coordinate,'native') || strcmpi(output_coordinate,'ACPC') || strcmpi(output_coordinate,'TAL') || strcmpi(output_coordinate,'MNI')
      fprintf('applying spatial transformation from FreeSurfer native MGZ to BrainVoyager VMR space...');
      trf=BVQXfile(fullfile(pwd,VMR_coregistration_file));
      vmr{hh}=mgzvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      vmr{hh}.ReferenceSpace=1;
      trf.ClearObject(); clear trf;
      fprintf('done.\n');
    end

    % transformation using the result of ACPC alignment
    if strcmpi(output_coordinate,'ACPC') || strcmpi(output_coordinate,'TAL') || strcmpi(output_coordinate,'MNI')
      fprintf('applying ACPC transformation...');
      trf=BVQXfile(fullfile(pwd,ACPC_file));
      tmpvmr=vmr{hh}.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      vmr{hh}.ClearObject();
      vmr{hh}=tmpvmr.CopyObject();
      vmr{hh}.ReferenceSpace=2; % ACPC coordinate
      tmpvmr.ClearObject(); clear tmpvmr;
      trf.ClearObject(); clear trf;
      fprintf('done.\n');
    end

    % transformation using the result of TAL conversion
    warning off; %#ok
    if strcmpi(output_coordinate,'TAL')
      fprintf('applying TAL transformation...');
      tal=BVQXfile(fullfile(pwd,TAL_file));

      % handling a new version of a TAL file generated by BrainVoayer 20 and above.
      fid=fopen(fullfile(pwd,TAL_file),'r');
      if fid==-1, error('can not open TAL_file. check the input variable.'); end
      curline=fgetl(fid);
      fclose(fid);
      if ~isempty(strfind(curline,'FileVersion'))
      %if isstructmember(tal,'REMAININGLINES') && length(tal.REMAININGLINES)>3
        newtal=BVQXfile('new:tal');
        idx=strfind(tal.REMAININGLINES{1},'"');
        newtal.Project=char(tal.REMAININGLINES{1}(idx(1)+1:end-1));
        val=textscan(tal.REMAININGLINES{2},'%s'); newtal.AC=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{3},'%s'); newtal.PC=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{4},'%s'); newtal.AP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{5},'%s'); newtal.PP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{6},'%s'); newtal.SP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{7},'%s'); newtal.IP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{8},'%s'); newtal.RP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        val=textscan(tal.REMAININGLINES{9},'%s'); newtal.LP=[str2num(val{1}{2}),str2num(val{1}{3}),str2num(val{1}{4})]; %#ok
        %newtal.REMAININGLINES(1,:)=tal.REMAININGLINES{10};
        %newtal.REMAININGLINES(2,:)=tal.REMAININGLINES{11};
        newtal.RunTimeVars=tal.RunTimeVars;
        tal.ClearObject(); clear tal;
        tal=newtal.CopyObject();
        newtal.ClearObject(); clear newtal;
      end
      tmpvmr=vmr{hh}.Talairach('linear',tal);
      vmr{hh}.ClearObject();
      vmr{hh}=tmpvmr.CopyObject();
      vmr{hh}.ReferenceSpace=3; % TAL coordinate
      vmr{hh}.VoxResVerified=1;
      tmpvmr.ClearObject(); clear tmpvmr;
      tal.ClearObject(); clear tal;
      fprintf('done.\n');
    end
    warning on; %#ok

    % transformation using the result of MNI transformation (BrainVoyager 20.2 or above is required)
    if strcmpi(output_coordinate,'MNI')
      fprintf('applying MNI transformation...');
      trf=BVQXfile(fullfile(pwd,MNI_file));
      tmpvmr=vmr{hh}.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      vmr{hh}.ClearObject();
      vmr{hh}=tmpvmr.CopyObject();
      vmr{hh}.ReferenceSpace=4; % MNI coordinate
      tmpvmr.ClearObject(); clear tmpvmr;
      trf.ClearObject(); clear trf;
      fprintf('done.\n');
    end

  else % if ~strcmpi(output_coordinate,'native')
    vmr{hh}=mgzvmr.CopyObject();
    vmr{hh}.ReferenceSpace=1; % native coordinate
  end

  % prepare for the surface reconstruction
  vmr{hh}.VMRData(vmr{hh}.VMRData>0)=240; % remap the value to avoid the interpolation effect.
  vmr{hh}=vmr{hh}.PrepareForReco(240,235);

  % clean up
  mgzvmr.ClearObject(); clear mgzvmr;
end % for hh=1:1:length(hemis)

% saving
fprintf('saving the converted VMR...');
for hh=1:1:length(hemis)
  vmr{hh}.SaveAs(fullfile(mgzpath,[mgzfname,mgzext,strrep(sprintf('_%s_%s_BL2.vmr',hemis{hh},output_coordinate),'_native','')]));
end
fprintf('done.\n');

% clean up
if ~nargout
  for hh=1:1:length(hemis)
    vmr{hh}.ClearObject();
  end
  clear vmr;
end

return

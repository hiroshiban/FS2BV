function vmr=ConvertFreeSurferMGZ2VMR(mgz_file,VMR_coregistration_file,ACPC_file,TAL_file,MNI_file,asis_conversion_flg,flip_flg,vmr_framingcube)

% Reads a FreeSurfer MGZ file and converts it to a BrainVoyager VMR file.
% function vmr=ConvertFreeSurferMGZ2MVR(mgz_file,:VMR_coregistration_file,:ACPC_file,:TAL_file,:MNI_file,:asis_conversion_flg,:flip_flg,:vmr_framingcube)
% (: is optional)
%
% This function reads a FreeSurfer MGZ file, such as aparc.a2009s+aseg.mgz (in ~/$FREESURFER_SUBJ_DIR/(subj)/mri/),
% and converts it to a BrainVoyager VMR file. If you just want to convert FreeSurfer MGZ to BrainVoyager VMR format,
% please just set mgz_file and run this function without any other input variables. If you want to generate the VMR
% files in the BrainVoyager-native subject anatomy space, please additionally set the 'VMR_coregistration_file'
% variable. If you want to create VMRs in ACPC, TAL, or MNI space, please further set 'ACPC_file' or/and 'TAL/MNI_file'
% respectively. To create TAL VMRs, you don't need to set MNI_file, while you don't need to set TAL_file (so please set
% TAL_file='';) to create MNI VMRs.
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
% >> mgz_file='./freesurfer/subjects/DC/aparc.a2009s+aseg.mgz';
% >> VMR_coregistration_file='./DC/3d/wm.seg-TO-DC_10_dc17_006.3d_final.trf';
% >> ACPC_file='./DC/3d/DC_10_dc17_006.3d_final_ACPC.trf';
% >> TAL_file='./DC/3d/DC_10_dc17_006.3d_final_TAL.tal';
% >> voi=ConvertFreeSurferMGZ2VMR(mgz_file,VMR_coregistration_file,ACPC_file,TAL_file,'',0);
%
% [input]
% mgz_file    : A FreeSurfer MGZ file to be converted to VMR specified with a relative path format,
%               in which the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC/aparc.a2009s+aseg.mgz';
% VMR_coregistration_file : (optional) a transformation matrix file that is generated after coregistering
%               (fs_subj_dir)/mri/wm.seg.mgz(vmr) to the corresponding subject's BrainVoyager VMR
%               file in the native space (generally *_final.vmr in BVQX_hbtools).
%               for the details of coregistration, please see MaskVMRbyFreeSurferSegmentation.m
%               or MaskVMRbyFreeSurferSegmentation_ribbon.m
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
% asis_conversion_flg : (optional) whether converting MGZ to VMR without converting the volume value.
%               1 by default.
%               If you are going to convert, for instance, an auto segmented result, e.g. aseg.mgz, to VMR,
%               it is required to adjust the values (IDs, specifying segmented regions/tissues, for details
%               please see FreeSurferColorLUT.txt) assined to voxels in the volume data since a MGZ stores
%               the volume values (IDs) in a float format (0-65535) while BrainVoyager VMR stores them as an
%               int8 (0-255) format. Therefore, a careful linear conversion from 0-65535 to 0-255 is required.
%               To do this, please set this input variable to 0. Otherwise, this function tries to convert
%               the values simply using out_img=(in_img-min(in_img))./max(in_img)*255;
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
% The generated VMR file is saved in the same directory with the input MGZ file with an extention, '.vmr'.
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
% Created    : "2017-08-01 10:19:14 ban"
% Last Update: "2021-06-16 10:02:41 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(VMR_coregistration_file), VMR_coregistration_file=''; end
if nargin<3 || isempty(ACPC_file), ACPC_file=''; end
if nargin<4 || isempty(TAL_file), TAL_file=''; end
if nargin<5 || isempty(MNI_file), MNI_file=''; end
if nargin<6 || isempty(asis_conversion_flg), asis_conversion_flg=1; end
if nargin<7 || isempty(flip_flg), flip_flg=1; end
if nargin<8 || isempty(vmr_framingcube), vmr_framingcube=256; end

% check the input files
if ~exist(fullfile(pwd,mgz_file),'file'), error('mgz_file not found. check the input variable.'); end

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
fprintf('Target FreeSurfer MGZ : %s\n',mgz_file);
fprintf('Output coordinate     : %s\n',output_coordinate);
fprintf('\n');

%% processing

% loading the input MGZ file
[mgzpath,mgzfname,mgzext]=fileparts(fullfile(pwd,mgz_file));
mgz=MRIread(fullfile(pwd,mgz_file));
%if uint8(mgz.volres(1))~=1
%  error('currently only 1mm cubic segmentation file is accepted. check the input FreeSurfer segmentation file.');
%end

% storing the original FreeSurfer segmentation result as a reference
try
  if ~flip_flg
    mgzvol=flip(flip(permute(mgz.vol,[3,1,2]),1),3);
  else
    mgzvol=flip(permute(mgz.vol,[3,1,2]),1);
  end
catch
  if ~flip_flg
    mgzvol=flipdim(flipdim(permute(mgz.vol,[3,1,2]),1),3);
  else
    mgzvol=flipdim(permute(mgz.vol,[3,1,2]),1); %#ok
  end
end

% volume value (region/tissue indices) conversions
fprintf('adjusting image intensities...');
if asis_conversion_flg
  % simple adjustment of the volume intensities from 0-65535 to 0-255
  %new_mgzvol=uint8(reshape((mgzvol(:)-min(mgzvol(:)))./(max(mgzvol(:))-min(mgzvol(:)))*255,size(mgzvol)));
  new_mgzvol=uint8(reshape((mgzvol(:)-min(mgzvol(:)))./(max(mgzvol(:))-min(mgzvol(:)))*225,size(mgzvol))); % as VMR values above 230 are reserved as some masks
else
  % special treatment to preserve region/tissue segmentation results
  %
  % BrainVoyager VMR stores pixel(voxel) gray-scale intensiteis as an 8 bit (0-255) format, while the FreeSurfer
  % MGZ stores in a float format (0-65526, allowing floating values). Therefore, we need to adjust the tissue IDs
  % in the MGZ file so that it ranges 0-255. Furthermore, in BrainVoyager VMR handling, 240 and the other values
  % are reserved for specific tissues such as white matter etc. I will adjust the volume IDs from 1 to 230.

  % get the number of regions/tissues (different values (IDs) in the target *.mgz file) included in the volume
  seg_idx=setdiff(unique(mgzvol),0);
  if numel(seg_idx)>230
    warning(['Too many regions/tissues are segmented and stored in the input MGZ file.',...
             'Some values may be rounded off in the next step. please be careful']);
  end

  new_seg_idx=round(linspace(1,230,numel(seg_idx)));
  new_mgzvol=zeros(size(mgzvol));
  for vv=1:1:numel(seg_idx)
    new_mgzvol(mgzvol==seg_idx(vv))=new_seg_idx(vv);
  end
end
fprintf('done.\n');

% generating an output vmr file
mgzvmr=BVQXfile('new:vmr');
mgzvmr.CoordinateSystem=2;
mgzvmr.ReferenceSpace=1; % native coordinate
mgzvmr.VoxResVerified=0;
mgzvmr.VoxResInTalairach=1;
mgzvmr.FramingCube=vmr_framingcube;
mgzvmr.DimX=mgz.height;
mgzvmr.DimY=mgz.width;
mgzvmr.DimZ=mgz.depth;
mgzvmr.NRows=mgzvmr.DimX;
mgzvmr.FOVRows=mgzvmr.DimX;
mgzvmr.FOVCols=mgzvmr.DimY;
mgzvmr.VMRData=uint8(new_mgzvol);

if ~strcmpi(output_coordinate,'native') || ~isempty(VMR_coregistration_file)

  % transformation using the result of VMR (wm.seg.vmr) to VMR (*_final.vmr) coregistration
  if ~isempty(VMR_coregistration_file)
  %if strcmpi(output_coordinate,'native') || strcmpi(output_coordinate,'ACPC') || strcmpi(output_coordinate,'TAL') || strcmpi(output_coordinate,'MNI')
    fprintf('applying spatial transformation from FreeSurfer native MGZ to BrainVoyager VMR space...');
    trf=BVQXfile(fullfile(pwd,VMR_coregistration_file));
    vmr=mgzvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
    vmr.ReferenceSpace=1;
    trf.ClearObject(); clear trf;
    fprintf('done.\n');
  end

  % transformation using the result of ACPC alignment
  if strcmpi(output_coordinate,'ACPC') || strcmpi(output_coordinate,'TAL') || strcmpi(output_coordinate,'MNI')
    fprintf('applying ACPC transformation...');
    trf=BVQXfile(fullfile(pwd,ACPC_file));
    tmpvmr=vmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
    vmr.ClearObject(); clear vmr;
    vmr=tmpvmr.CopyObject();
    vmr.ReferenceSpace=2; % ACPC coordinate
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
    tmpvmr=vmr.Talairach('linear',tal);
    vmr.ClearObject(); clear vmr;
    vmr=tmpvmr.CopyObject();
    vmr.ReferenceSpace=3; % TAL coordinate
    vmr.VoxResVerified=1;
    tmpvmr.ClearObject(); clear tmpvmr;
    tal.ClearObject(); clear tal;
    fprintf('done.\n');
  end
  warning on; %#ok

  % transformation using the result of MNI transformation (BrainVoyager 20.2 or above is required)
  if strcmpi(output_coordinate,'MNI')
    fprintf('applying MNI transformation...');
    trf=BVQXfile(fullfile(pwd,MNI_file));
    tmpvmr=vmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
    vmr.ClearObject(); clear vmr;
    vmr=tmpvmr.CopyObject();
    vmr.ReferenceSpace=4; % MNI coordinate
    tmpvmr.ClearObject(); clear tmpvmr;
    trf.ClearObject(); clear trf;
    fprintf('done.\n');
  end

else % if ~strcmpi(output_coordinate,'native')
  vmr=mgzvmr.CopyObject();
  vmr.ReferenceSpace=1; % native coordinate
end

% saving
fprintf('saving the converted VMR...');
vmr.SaveAs(fullfile(mgzpath,[mgzfname,mgzext,strrep(sprintf('_%s.vmr',output_coordinate),'_native','')]));
fprintf('done.\n');

% clean up
mgzvmr.ClearObject(); clear mgzvmr;
if ~nargout
  vmr.ClearObject(); clear vmr;
end

return

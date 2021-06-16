function voi=ConvertFreeSurferParcellation2BVvoi(fs_subj_dir,VMR_coregistration_file,ACPC_file,TAL_file,MNI_file,flilp_flg,vmr_framingcube)

% Reads FreeSurfer auto parcellation files and converts them to BrainVoyager VOI files.
% function voi=ConvertFreeSurferParcellation2BVvoi(fs_subj_dir,:VMR_coregistration_file,:ACPC_file,:TAL_file,:MNI_file,:flip_flg,:vmr_framingcube)
% (: is optional)
%
% This function reads FreeSurfer auto parcellation files, aparc.a2009s+aseg.mgz, aparc.DKTatlas+aseg.mgz,
% aparc+aseg.mgz, and aseg.mgz (those files are in ~/fs_subj_dir/mri/), and converts them to BrainVoyager
% VOI files. If you want to simply convert auto segmentation files to VOIs without any transformation,
% please just set fs_suj_dir alone. If you want to generate the VOI files in the BrainVoyager-native
% subject anatomy space, please additionally set the 'VMR_coregistration_file' variable. If you want to
% create VOIs in ACPC, TAL, or MNI space, please further set 'ACPC_file' or/and 'TAL/MNI_file' respectively.
% To create TAL VOIs, you don't need to set MNI_file, while you don't need to set TAL_file (so please set
% TAL_file='';) to create MNI VOIs.
%
% For FreeSurfer details, please see
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki
% https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
% https://surfer.nmr.mgh.harvard.edu/fswiki/SubcorticalSegmentation
% https://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011
%
% !!!IMPORTANT NOTE!!!
% It may take so much computation time to convert FreeSurfer ROIs into VOI in TAL space using this function.
% Therefore it would be recommended to create VOI files in ACPC space and convert them into TAL space
% using the transformation option implemented in BrainVoyager's Volume-Of-Interest Analysis GUI window.
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
% >> fs_subj_dir='./freesurfer/subjects/DC';
% >> VMR_coregistration_file='./DC/3d/wm.seg-TO-DC_10_dc17_006.3d_final.trf';
% >> ACPC_file='./DC/3d/DC_10_dc17_006.3d_final_ACPC.trf';
% >> TAL_file='./DC/3d/DC_10_dc17_006.3d_final_TAL.tal';
% >> voi=ConvertFreeSurferParcellation2BVvoi(fs_subj_dir,VMR_coregistration_file,ACPC_file,TAL_file);
%
% [input]
% fs_subj_dir : FreeSurfer subject directory with a relative path format, in which
%               the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC
% VMR_coregistration_file : (optional) a transformation matrix file that is generated after
%               coregistering (fs_subj_dir)/mri/wm.seg.mgz(vmr) to the corresponding subject's
%               BrainVoyager VMR file in the native space (generally *_final.vmr in BVQX_hbtools).
%               for the details of coregistration, please see MaskVMRbyFreeSurferSegmentation.m
%               The file should be specified as a relative path format.
%               If this variable is not set, the auto segmentation files are simply converted
%               from *.mgz to BrainVoyager VOI files without any transformation.
%               empty by default.
% ACPC_file   : (optional) BrainVoyager ACPC transformation file (*_ACPC.trf), specified
%               as a relative path format. If not specified, VOIs are generated in the native
%               space. empty by default.
% TAL_file    : (optional) BrainVoyager TAL transformation file (*_TAL.tal), specified
%               as a relative path format. If not specified, VOIs are generated in the native
%               or ACPC space. To generate a VOI file in TAL space, please also specify ACPC_file.
%               empty by default.
% MNI_file    : (optional) BrainVoyager MNI transformation file (*_TO_MNI_a12.trf), specified
%               as a relative path format. If not specified, VOIs are generated in the native,
%               ACPC, or TAL space. To generate a VOI file in MNI space, please also specify
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
% voi         : BrainVoyager VOI object, a cell structure
% The generated VOI files are saved in the same directory with the input auto parcellation files
% with an extention '.voi'.
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
if nargin<6 || isempty(flip_flg), flip_flg=1; end
if nargin<7 || isempty(vmr_framingcube), vmr_framingcube=256; end

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
if ~isempty(ACPC_file), output_coordinate='ACPC'; end
if ~isempty(TAL_file), output_coordinate='TAL'; end
if ~isempty(MNI_file), output_coordinate='MNI'; end

%% check the freesurfer directory and files
if ~exist(fullfile(pwd,fs_subj_dir),'dir')
  error('fs_subj_dir not found. check the input variable.');
end

% set files to be processed
% note: aseg is subcortical segmentation results, and the others are cortical gray matter segmentations
fsfiles{1}=fullfile(pwd,fs_subj_dir,'mri','aparc.a2009s+aseg.mgz');
fsfiles{2}=fullfile(pwd,fs_subj_dir,'mri','aparc.DKTatlas+aseg.mgz');
fsfiles{3}=fullfile(pwd,fs_subj_dir,'mri','aparc+aseg.mgz');
fsfiles{4}=fullfile(pwd,fs_subj_dir,'mri','aseg.mgz');

for ii=1:1:length(fsfiles)
  if ~exist(fsfiles{ii},'file')
    %[dummy,fsfname,fsext]=fileparts(fsfiles{ii});
    %warning('%s%s not found. omitting it from the processing',fsfname,fsext);
    fsfiles{ii}='';
  end
end

% set the path to the tools for handling FreeSurfer files on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

% display message
fprintf('Target FreeSurfer subject directory : %s\n',fs_subj_dir);
fprintf('Output coordinate                   : %s\n',output_coordinate);
fprintf('\n');

% loading FreeSurfer ColorLUT definitions
[lut_dict,lutIDs]=readFreeSurferColorLUT();

%% processing
roival=150;
voi=cell(length(fsfiles),1);
for ii=1:1:length(fsfiles)

  if isempty(fsfiles{ii}), continue; end

  % display message
  [fspath,fsfname,fsext]=fileparts(fsfiles{ii});
  fprintf('converting FreeSurfer segmentation file to BrainVoayger VMR: %s%s...',fsfname,fsext);

  % converting *.mgz to *.vmr
  fsfile=MRIread(fsfiles{ii});
  %if ~isempty(find(uint8(fsfile.volres(1))~=1,1))
  %  error('currently only 1mm cubic segmentation file is accepted. check the input FreeSurfer segmentation file.');
  %end
  if numel(unique(fsfile.volsize))~=1
    error('currently only cubic volume (size_x = size_y = size_z) datasets is accepted. check the input FreeSurfer segmentation file.');
  end

  fprintf('done.\n');

  % storing the original FreeSurfer segmentation result as a reference
  try
    if ~flip_flg
      segresults=flip(flip(permute(fsfile.vol,[3,1,2]),1),3);
    else
      segresults=flip(permute(fsfile.vol,[3,1,2]),1);
    end
  catch
    if ~flip_flg
      segresults=flipdim(flipdim(permute(fsfile.vol,[3,1,2]),1),3);
    else
      segresults=flipdim(permute(fsfile.vol,[3,1,2]),1); %#ok
    end
  end

  % IMPORTANT NOTEs
  % from here, we have to run the VOI coordinate conversion one by one with adjusting voxel IDs.
  % otherwise, the default voxel values (double precision, around 0-12000) are rounded off to
  % uint8 (0-255) in saving, due to BrainVoyager VMR format limitations. Then, the VOI extractions will fail.
  fprintf('converting FreeSurfer ROIs to BrainVoyager VOIs in %s space...\n',output_coordinate);
  voicounter=0;

  % extract the LUT IDs stored in the current FreeSurfer segmentation volume and get the
  % corresponding LUT data structure from lut_dict. then we can save computation time as
  % it is enough to process only those ROIs in the following steps without checking the
  % whole LUT contents
  cur_lutIDs=unique(segresults(:));
  [dummy,IA]=intersect(lutIDs,cur_lutIDs);
  lut=lut_dict(IA);

  % initializing VOI file
  voi{ii}=BVQXfile('new:voi');
  voi{ii}.FileVersion=4; % for compatibility with the latest BrainVoyager 20.6, I fixe the VOI file vertion as 4.
  voi{ii}.OriginalVMRFramingCubeDim=vmr_framingcube;

  % generating a vmr file as a workspace
  fsvmr=BVQXfile('new:vmr');
  fsvmr.FramingCube=vmr_framingcube;
  fsvmr.DimX=fsfile.height;
  fsvmr.DimY=fsfile.width;
  fsvmr.DimZ=fsfile.depth;
  fsvmr.NRows=fsvmr.DimX;
  fsvmr.FOVRows=fsvmr.DimX;
  fsvmr.FOVCols=fsvmr.DimY;
  fsvmr.VMRData=zeros(size(segresults));

  % processing
  for vv=1:1:length(lut)

    if strcmpi(lut{vv}.name,'Unknown'), continue; end

    % display message
    voicounter=voicounter+1;
    fprintf('processing #ROI %04d (ID: %5d): %s...',voicounter,lut{vv}.ID,lut{vv}.name);

    % extracting only the target ROI to avoid rounding off problem described above and
    % update fsvmr for the current ROI, with setting a specific volumes defined above.
    fsvmr.VMRData(:)=0; % clean up with saving memory allocation time
    if ~strcmpi(output_coordinate,'native')
      fsvmr.VMRData(segresults==lut{vv}.ID)=roival;
    else
      voi_idx=find(segresults==lut{vv}.ID);
      fsvmr.VMRData(voi_idx)=roival;
    end

    % transformation using the result of VMR (wm.seg.vmr) to VMR (*_final.vmr) coregistration
    if ~isempty(VMR_coregistration_file)
      trf=BVQXfile(fullfile(pwd,VMR_coregistration_file));
      newvmr=fsvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      trf.ClearObject(); clear trf;
    elseif ~strcmpi(output_coordinate,'native')
      newvmr=fsvmr.CopyObject();
    end

    % transformation using the result of ACPC alignment
    if strcmpi(output_coordinate,'ACPC') || strcmpi(output_coordinate,'TAL') || strcmpi(output_coordinate,'MNI')
      trf=BVQXfile(fullfile(pwd,ACPC_file));
      tmpvmr=newvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      newvmr.ClearObject(); clear newvmr;
      newvmr=tmpvmr.CopyObject();
      tmpvmr.ClearObject(); clear tmpvmr;
      trf.ClearObject(); clear trf;
    end

    % transformation using the result of TAL registration matrix
    warning off; %#ok
    if strcmpi(output_coordinate,'TAL')
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
      tmpvmr=newvmr.Talairach('linear',tal);
      newvmr.ClearObject(); clear newvmr;
      newvmr=tmpvmr.CopyObject();
      tmpvmr.ClearObject(); clear tmpvmr;
      tal.ClearObject(); clear tal;
    end
    warning on; %#ok

    % transformation using the result of MNI registration matrix (BrainVoyager 20.2 or above is required)
    if strcmpi(output_coordinate,'MNI')
      trf=BVQXfile(fullfile(pwd,MNI_file));
      tmpvmr=newvmr.ApplyTRF(trf,struct('asdouble',true,'inverse',false,'method','nearest'));
      newvmr.ClearObject(); clear newvmr;
      newvmr=tmpvmr.CopyObject();
      tmpvmr.ClearObject(); clear tmpvmr;
      trf.ClearObject(); clear trf;
    end

    % get the voxel indices which exceeds the 'thres'
    if ~strcmpi(output_coordinate,'native')
      voi_idx=find(newvmr.VMRData>0);%find( roival-10<=newvmr.VMRData & newvmr.VMRData<=roival+10 );
      voi_voxels=zeros(numel(voi_idx),3);
      [voi_voxels(:,1),voi_voxels(:,2),voi_voxels(:,3)]=ind2sub(size(newvmr.VMRData),voi_idx);
    else
      voi_voxels=zeros(numel(voi_idx),3);
      [voi_voxels(:,1),voi_voxels(:,2),voi_voxels(:,3)]=ind2sub(size(fsvmr.VMRData),voi_idx);
    end

    % here -1 is required probably as FreeSurfer coordinate starts from 1 while BrainVoayger starts from 0 (have to check).
    % or probably as volume starting points are different between FreeSurfer and BrainVoyager.
    voi_voxels=voi_voxels-1;

    % set the new VOI
    voi{ii}.NrOfVOIs=voi{ii}.NrOfVOIs+1;
    voi{ii}.VOI(voi{ii}.NrOfVOIs).Name=lut{vv}.name;
    voi{ii}.VOI(voi{ii}.NrOfVOIs).Color=lut{vv}.RGB;
    if strcmpi(output_coordinate,'native') || strcmpi(output_coordinate,'ACPC')
      voi{ii}.VOI(voi{ii}.NrOfVOIs).NrOfVoxels=size(voi_voxels,1);
      voi{ii}.VOI(voi{ii}.NrOfVOIs).Voxels=voi_voxels;
    else % strcmpi(output_coordinate,'TAL')
      voi{ii}.VOI(voi{ii}.NrOfVOIs).NrOfVoxels=size(voi_voxels,1);
      voi{ii}.VOI(voi{ii}.NrOfVOIs).Voxels=round(fsfile.volsize(1)/2)-voi_voxels(:,[3,1,2]); % not sure but special coordinate transformations are required for TAL VOIs
    end

    fprintf('done.\n');

  end % for vv=1:1:length(lut)

  % clean up
  fsvmr.ClearObject(); clear fsvmr;

  % saving the result
  fprintf('saving the generated VOI file...');
  if strcmpi(output_coordinate,'native')
    voi{ii}.ReferenceSpace='BV';
    %voi_prefix='_native';
    voi_prefix='';
  elseif strcmpi(output_coordinate,'ACPC')
    voi{ii}.ReferenceSpace='ACPC';
    voi_prefix='_ACPC';
  elseif strcmpi(output_coordinate,'TAL')
    voi{ii}.ReferenceSpace='TAL';
    voi_prefix='_TAL';
  else
    voi_prefix='';
  end

  % saving the generated VOI file
  %voi{ii}.SaveAs(fullfile(fspath,[fsfname,fsext,voi_prefix,'.voi']));

  % here I manually write the VOI file to save the file size (if we use voi.SaveAs() function, redundant white spaces are inserted in the file)
  fid=fopen(fullfile(fspath,[fsfname,fsext,voi_prefix,'.voi']),'W');
  if fid==-1, error('%s can not be open. check the file, run fclose all, and try again.',[fsfname,fsext,voi_prefix,'.voi']); end

  % write the VOI details into the file
  fprintf(fid,'\nFileVersion:                %d\n',voi{ii}.FileVersion);

  if voi{ii}.FileVersion<4
    fprintf(fid,'\nCoordsType:                 %s\n',voi{ii}.ReferenceSpace);
    fprintf(fid,'\nSubjectVOINamingConvention: %s\n',voi{ii}.SubjectVOINamingConvention);
    fprintf(fid,'\n');
  else
    fprintf(fid,'\nReferenceSpace:             %s\n',voi{ii}.ReferenceSpace);
    fprintf(fid,'\n');
    fprintf(fid,'OriginalVMRResolutionX:     %d\n',voi{ii}.OriginalVMRResolutionX);
    fprintf(fid,'OriginalVMRResolutionY:     %d\n',voi{ii}.OriginalVMRResolutionY);
    fprintf(fid,'OriginalVMRResolutionZ:     %d\n',voi{ii}.OriginalVMRResolutionZ);
    fprintf(fid,'OriginalVMROffsetX:         %d\n',voi{ii}.OriginalVMROffsetX);
    fprintf(fid,'OriginalVMROffsetY:         %d\n',voi{ii}.OriginalVMROffsetY);
    fprintf(fid,'OriginalVMROffsetZ:         %d\n',voi{ii}.OriginalVMROffsetZ);
    fprintf(fid,'OriginalVMRFramingCubeDim:  %d\n',voi{ii}.OriginalVMRFramingCubeDim);
    fprintf(fid,'\nLeftRightConvention:        %d\n',voi{ii}.Convention);
    fprintf(fid,'\nSubjectVOINamingConvention: %s\n',voi{ii}.SubjectVOINamingConvention);
  end

  fprintf(fid,'\nNrOfVOIs:                   %d\n',voi{ii}.NrOfVOIs);
  fprintf(fid,'\n');

  for vv=1:1:voi{ii}.NrOfVOIs
    fprintf(fid,'NameOfVOI:                  %s\n',voi{ii}.VOI(vv).Name);
    fprintf(fid,'ColorOfVOI:                 %d %d %d\n',(voi{ii}.VOI(vv).Color)');
    fprintf(fid,'\n');
    fprintf(fid,'NrOfVoxels:                 %d\n',voi{ii}.VOI(vv).NrOfVoxels);
    voxelfmt='%d %d %d\n'; % * HERE, REMOVING THE REDUNDANT WHITE SPACES *
    fprintf(fid,voxelfmt,(voi{ii}.VOI(vv).Voxels)');
    fprintf(fid,'\n\n');
  end
  fprintf(fid,'NrOfVOIVTCs:                %d',voi{ii}.NrOfVTCs);

  % fclose the VOI file
  fclose(fid);

  fprintf('done.\n');
  fprintf('\n');

  % clean up
  if ~nargout
    if ~isempty(voi{ii})
      voi{ii}.ClearObject();
    end
  end

end % for ii=1:1:length(fsfiles)

% remove the path to the tools for handling FreeSurfer files on MATLAB
rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

return

function nii=ExtractFSLroi(XML_dir,ROI_label_idx,ROI_name,prefix_xml,save_flg)

% Reads ROI XML file(s) defined by FSL software (NIFTI file) and extracts the target ROI as a separate NII file.
% function nii=ExtractFSLroi(XML_dir,label_idx,ofname,:prefix_xml,:save_flg)
% (: is optional)
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
% To find the VOIs in which specific XYZ coordinates in TAL/MNI space, please use the function below,
% GetAreaNameFromAtlasVOI               : Returns area candidates, in which the input XYZ coordinate(s)
%                                        is(are) belonging to, based on the pre-defined VOI atlases.
%
% Some *.nii/*.nii.gz/*.img files are provided as ROIs 2mm cubic volume space.
% If you want to convert them to 1mm cubic space, that is compatible with BrainVoayger default anatomical
% volume space, please run the command below
% >> reslice_nii('AAL.nii','AAL_1mm.nii',1,'','',2);
% before using the functions above. 'reslice_nii' is in BVQX_hbtools/nifti_tools.
%
% [input]
% XML_dir       : directory that contains *.xml files describing ROI(s) details defined in the *.nii.gz files
%                 a relative path format. the relative path should be set so that the location where this function
%                 is called is the origin. e.g. XML_dir='../../ROI_nii'
% ROI_label_idx : ROI label index defined in XML.
%                 e.g. if the target XML is Talairach.xml and the target ROI is 'Pulvinar', please find the line below
%                 in the target XML file.
%                 .....
%                 <label index="622" x="52" y="49" z="39">Left Cerebrum.Sub-lobar.Thalamus.Gray Matter.Pulvinar</label>
%                 .....
%                 then, please set label_idx=622;
%                 multiple values are acceptable. e.g. label_idx=[622,623];
% ROI_name      : output file name prefix, e.g. ofname='Pulvinar'; you do not need to include '_' and a file extension.
% prefix_xml    : (optional) file prefix to specify the target XML from multiple files.
%                 e.g. prefix_xml='Talairach', empty by default.
% save_flg      : (optional) whether saving the extracted ROI as NII file. [0|1]
%                 1 by default. when the target ROI nii_file is 'ABC.nii.gz', then the ROI will be saved as
%                 'ABC_{ROI_name}.nii.gz' in the directory where this function is called. when the target is 'ABC.nii',
%                 the ROI will be saved as non-compressed NII file as 'ABC_{ROI_name}.nii'.
%                 In the saved file, the target ROI region will be marked as 1 while the background is 0.
%
% [output]
% nii        : NIFTI file header and data, a matlab structure
%
% [example]
% >> %% some constants
% >> XML_dir='./FSL_atlases';
% >>
% >> % in the target XML file, Pulvinar indeces are
% >> % <label index="622" x="52" y="49" z="39">Left Cerebrum.Sub-lobar.Thalamus.Gray Matter.Pulvinar</label>
% >> % <label index="623" x="35" y="49" z="39">Right Cerebrum.Sub-lobar.Thalamus.Gray Matter.Pulvinar</label>
% >> ROI_label_idx{1}=[622,623];
% >> ROI_name{1}='Pulvinar';
% >>
% >> % in the target XML file, LGN indeces are
% >> % <label index="497" x="55" y="50" z="34">Left Cerebrum.Sub-lobar.*.Gray Matter.Lateral Geniculum Body</label>
% >> % <label index="500" x="32" y="51" z="34">Right Cerebrum.Sub-lobar.*.Gray Matter.Lateral Geniculum Body</label>
% >> ROI_label_idx{2}=[497,500];
% >> ROI_name{2}='LGN';
% >>
% >> prefix_xml='Talairach';
% >> save_flg=1;
% >>
% >> %% processing
% >> for ii=1:1:length(ROI_label_idx), ExtractFSLroi(XML_dir,ROI_label_idx{ii},ROI_name{ii},prefix_xml,save_flg); end
%
% [note on how to set the 'prefix_*' variable]
% prefix_* can be set flexibly as below.
% 1. a string: setting an including prefix (string) alone
%    e.g. prefix_*='_TDTS6.0';
%         --> processes files whose names contain '_TDTS6.0'
% 2. a {1 x N} cell string: setting including prefix (string) arrays
%    e.g. prefix_*={'_TDTS6.0','_TSS5.0mm'};
%         --> processes files whose names contain '_TDTS6.0s' or '_TSS5.0mm'.
% 3. a {2 x N} cell string: setting including/excluding prefix (string) arrays
%    e.g. prefix_*={{'_TDTS6.0s','_TSS5.0mm'};{'THP'}};
%         --> processes files whose names contain '_TDTS6.0s'
%             or '_TSS5.0mm' but do not contain 'THP'.
%         prefix_*={'';{'_TDTS6.0s'}};
%         --> processes files whose names do not contain '_TDTS6.0s'.
%         prefix_*={'_TSS5.0mm';''};
%         --> processes files whose names contain '_TSS5.0mm'.
%
% [reference]
% FSL by FMRIB, Oxford
% ref: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
%
% [dependency]
% 1. BVQXtools v0.8d
% ref: http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools/232-getting-started.html
%
% 2. BVQX_hbtools by H.Ban
%
% 3. nii_tools (all required functions are already included as subfunctions in this file)
% Developed by Jimmy Shen (pls@rotman-baycrest.on.ca), Modified by H.Ban for this function
%
%
% Created    : "2015-12-26 11:03:16 ban"
% Last Update: "2018-06-09 17:19:25 ban"

% check the input variables
if nargin<3, help(mfilename()); return; end
if nargin<4 || isempty(prefix_xml), prefix_xml=''; end
if nargin<5 || isempty(save_flg), save_flg=1; end

% check the input nii directory
XML_dir=fullfile(pwd,XML_dir);
if ~exist(XML_dir,'dir'), error('can not find XML_dir. check input variables.'); end
fprintf('Target XML dir: %s\n',XML_dir);

% get XML files
xmlfiles=GetFiles(XML_dir,'*.xml',prefix_xml);
if isempty(xmlfiles), error('can not find any XML file in XML_dir. check the input variables.'); end

% processing
for ii=1:1:length(xmlfiles)
  [xmlpath,xmlfname,xmlext]=fileparts(xmlfiles{ii});
  fprintf('\n');
  fprintf('Processing: %s%s...\n',xmlfname,xmlext);

  % parse the target XML contents
  xml=xmlread(xmlfiles{ii});

  % the codes to get NII file prefix below are messy, just generated based on the current FSL data conventions

  % get the corresponding NII file prefix
  try
    nii_prefix=char(xml.getElementsByTagName('header').item(0)...
                    .getElementsByTagName('shortname').item(0)...
                    .getFirstChild().getData());
  catch
    nii_prefix='';
  end

  % get another NII file prefix to avoid error
  try
    nii_prefix_sub=char(xml.getElementsByTagName('header').item(0)...
                           .getElementsByTagName('images').item(0)...
                           .getElementsByTagName('imagefile').item(0)...
                           .getFirstChild().getData());
    [dummy,nii_prefix_sub]=fileparts(nii_prefix_sub); %#ok
    sub_idx=strfind(nii_prefix_sub,'-');
    nii_prefix_sub=nii_prefix_sub(1:sub_idx(2)-1);
  catch
    nii_prefix_sub='';
  end

  % get the data directory
  tgt_dir=fileparts(char(xml.getElementsByTagName('header').item(0)...
                            .getElementsByTagName('images').item(0)...
                            .getElementsByTagName('imagefile').item(0)...
                            .getFirstChild().getData()));

  % get the ROI data
  ROIdata=xml.getElementsByTagName('label');
  nROIs=ROIdata.getLength();
  ROIs=cell(1,1);
  roi_counter=0;
  for rr=1:1:nROIs
    tgt_roi=ROIdata.item(rr-1); % as item starts from 0
    tmp_idx=str2num(char(tgt_roi.getAttribute('index'))); %#ok

    if ~isempty(intersect(tmp_idx,ROI_label_idx))
      roi_counter=roi_counter+1;
      ROIs{roi_counter}.name=char(tgt_roi.getFirstChild.getData());
      if strcmp(unique(char(ROIdata.item(0).getFirstChild.getData())),'*.')
        % in some case, the first index (0) is like *.*.*.*.* which means the background
        % then we need to set index as it is.
        ROIs{roi_counter}.index=str2num(char(tgt_roi.getAttribute('index'))); %#ok
      else
        % otherwise we need to set the index with +1 as 0 corresponds to the background in the target NII file.
        ROIs{roi_counter}.index=str2num(char(tgt_roi.getAttribute('index')))+1; %#ok
      end
      ROIs{roi_counter}.x=str2num(char(tgt_roi.getAttribute('x'))); %#ok
      ROIs{roi_counter}.y=str2num(char(tgt_roi.getAttribute('y'))); %#ok
      ROIs{roi_counter}.z=str2num(char(tgt_roi.getAttribute('z'))); %#ok
      fprintf('  [%03d] % 20s (x,y,z)=(%d,%d,%d)\n',...
              ROIs{roi_counter}.index,ROIs{roi_counter}.name,...
              ROIs{roi_counter}.x,ROIs{roi_counter}.y,ROIs{roi_counter}.z);
    end
  end

  % check whether the all indeces are included before the latter ROI extraction procedures.
  tmp_index=[];
  for rr=1:1:length(ROIs), tmp_index=[tmp_index,ROIs{rr}.index]; end %#ok
  check_idx=setdiff(ROI_label_idx,tmp_index);
  if ~isempty(check_idx)
    for mm=1:1:numel(check_idx)
      fprintf('WARNING: ROI label %d not found. skipping...\n',check_idx(mm));
    end
  end

  fprintf('\n');

  % get all the corresponding ROI NII files.
  % NOTE: the codes below are messy but required to avoid the errors due to the differences of
  %       characters used in the database XML and the actual file ('_' and '-')
  niifiles=GetFiles(fullfile(xmlpath,tgt_dir),'*.nii*',nii_prefix);
  if isempty(niifiles)
    niifiles=GetFiles(fullfile(xmlpath,tgt_dir),'*.nii*',strrep(nii_prefix,'_','-'));
    if isempty(niifiles)
      niifiles=GetFiles(fullfile(xmlpath,tgt_dir),'*.nii*',strrep(nii_prefix,'-','_'));
      if isempty(niifiles)
        niifiles=GetFiles(fullfile(xmlpath,tgt_dir),'*.nii*',nii_prefix_sub);
        if isempty(niifiles)
          error('can not find any NII file in the target directory. check the input variables.');
        end
      end
    end
  end

  % extracting ROI data
  for mm=1:1:length(niifiles)
    [dumny,niifname,niiext]=fileparts(niifiles{mm}); %#ok
    niifile=[niifname,niiext];
    if ( length(niifile)>4 && strcmp(niifile(end-3:end),'.nii') ), niiext='.nii'; end
    if ( length(niifile)>7 && strcmp(niifile(end-6:end),'.nii.gz') ), niiext='.nii.gz'; end

    % skip if some irrelevant file is unexpectedly selected in the codes above.
    if ~strcmpi(niiext,'.nii') && ~strcmpi(niiext,'.nii.gz'), continue; end

    % load the nii file
    fprintf('extracting: %s%s --> %s_%s%s...',niifname,niiext,strrep(niifname,'.nii',''),ROI_name,niiext);
    if strcmpi(niiext,'.nii')
      nii=load_nii(niifiles{mm});
    elseif strcmpi(niiext,'.nii.gz')
      nii=load_niigz(niifiles{mm});
    else
      error('can not load nii file. acceptable format is one of ''*.nii'' and ''*.nii.gz''. check the input file');
    end

    % selecting only the first volume when nii.img is 4D data strucutre just for simplicity.
    % here I will select the first volume which generally contains the lower probability ROIs = more voxels.
    if numel(size(nii.img))==4, nii.img=squeeze(nii.img(:,:,:,1)); end

    % get the voxel indeces and set them to the BVQX VOI objects
    voi_idx=[];
    for rr=1:1:length(ROIs)
      if max(nii.img(:))==2*nROIs
        % in some case, ROI labels in XML files are like index 1: area 1, index 2: area 2, ..., while the actual indeces are
        % like index 1: area 1 in the left hemisphere, index 2: the same area 1 in the right hemisphere, index 3: area 2 in the left hemisphere...
        % To handle that indexing converntion, we need to set the ROI indeces so that one area covers 2*index-1 & 2*index values
        voi_idx=[voi_idx;find(nii.img==2*ROIs{rr}.index-1 | nii.img==2*ROIs{rr}.index)]; %#ok
      else
        % otherwise, we can set the ROI index as it is
        voi_idx=[voi_idx;find(nii.img==ROIs{rr}.index)]; %#ok
      end
    end

    % masking nii.img and extract ROI(s)
    nii.img=zeros(size(nii.img));
    nii.img(voi_idx)=1;
    fprintf('done.\n');

    % save the new ROI NII file
    if save_flg
      fprintf('saving...');
      if strcmpi(niiext,'.nii')
        save_nii(nii,fullfile(pwd,[niifname,'_',ROI_name,'.nii']));
      elseif strcmpi(niiext,'.nii.gz')
        save_niigz(nii,fullfile(pwd,[strrep(niifname,'.nii',''),'_',ROI_name,'.nii.gz']));
      end
      fprintf('done.\n');
    end

  end % for mm=1:1:length(niifiles)
end % for ii=1:1:length(xmlfiles)

return


%% subfunctions

%---------------------------------------------------------------------
%  scaling etc. Other transforms (any degree rotation, shears, etc.)
%  are not supported, because in those transforms, each voxel has
%  to be repositioned, interpolated, and whole image(s) will have
%  to be reconstructed. If an input data (nii) can not be handled,
%  The program will exit with an error message "Transform of this
%  NIFTI data is not supported by the program". After the transform,
%  nii will be in RAS orientation, i.e. X axis from Left to Right,
%  Y axis from Posterior to Anterior, and Z axis from Inferior to
%  Superior. The RAS orientation system sometimes is also referred
%  as right-hand coordinate system, or Neurologist preferred system.
%
%  Usage: [nii] = load_nii(filename, [img_idx], [old_RGB])
%
%  filename - NIFTI file name.
%
%  img_idx    - a numerical array of image indices. Only the specified
% images will be loaded. If there is no img_idx, all available
% images will be loaded.
%
%  old_RGB    - a boolean variable to tell difference from new RGB24 from old
%       RGB24. New RGB24 uses RGB triple sequentially for each voxel, like
%       [R1 G1 B1 R2 G2 B2 ...]. Analyze 6.0 developed by AnalyzeDirect uses
%       old RGB24, in a way like [R1 R2 ... G1 G2 ... B1 B2 ...] for each
%       slices. If the image that you view is garbled, try to set old_RGB
%       variable to 1 and try again, because it could be in old RGB24.
%
%  The number of images scans can be obtained from get_nii_frame, or
%  simply: hdr.dime.dim(5)
%
%  Returned values:
%
%  nii.hdr - struct with NIFTI header fields.
%  nii.filetype - Analyze format (0); NIFTI .hdr/.img (1); NIFTI .nii (2)
%  nii.fileprefix - NIFTI filename without extension.
%  nii.machine - machine string variable.
%  nii.img_idx - Indices of images to be loaded.
%  nii.img - 3D (or 4D) matrix of NIFTI data.
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%

function [nii] = load_nii(filename, img_idx, old_RGB)

if ~exist('filename','var'),
  error('Usage: [nii] = load_nii(filename, [img_idx], [old_RGB])');
end

if ~exist('img_idx','var'), img_idx = []; end
if ~exist('old_RGB','var'), old_RGB = 0; end

%  Read the dataset header
%
[nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

%  Read the dataset body
%
[nii.img,nii.hdr] = ...
    load_nii_img(nii.hdr,nii.filetype,nii.fileprefix,nii.machine,img_idx,old_RGB);

%  Perform some of sform/qform transform
%
nii = xform_nii(nii);

return         % load_nii


%---------------------------------------------------------------------
%  Load NIFTI dataset body after its header is loaded using load_nii_hdr.
%
%  Usage: [img,hdr] = ...
%       load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[old_RGB]);
%
%  Where: [hdr,filetype,fileprefix,machine] = load_nii_hdr(filename);
%
%  img_idx    - a numerical array of image indices. Only the specified images
% will be loaded. If there is no img_idx, all images will be loaded. The
%       number of images scans can be obtained from hdr.dime.dim(5)
%
%  old_RGB    - an optional boolean variable to handle special RGB data
%       sequence [R1 R2 ... G1 G2 ... B1 B2 ...] that is used only by
%       AnalyzeDirect (Analyze Software). Since both NIfTI and Analyze
%       file format use RGB triple [R1 G1 B1 R2 G2 B2 ...] sequentially
%       for each voxel, this variable is set to FALSE by default. If the
%       image that you displayed is garbled, try to set old_RGB variable
%       to TRUE (or 1) and load the image again, because it could be using
%       AnalyzeDirect RGB data sequence.
%
%  Returned values:
%
%  img - 3D (or 4D) matrix of NIFTI data.
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%

function [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,img_idx,old_RGB)

if ~exist('hdr','var') || ~exist('filetype','var') || ~exist('fileprefix','var') || ~exist('machine','var')
  error('Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx]);');
end

if ~exist('img_idx','var') || hdr.dime.dim(5)<1, img_idx = []; end
if ~exist('old_RGB','var'), old_RGB = 0; end

%  check img_idx
%
if ~isempty(img_idx) && ~isnumeric(img_idx)
  error('"img_idx" should be a numerical array.');
end

if length(unique(img_idx)) ~= length(img_idx)
  error('Duplicate image index in "img_idx"');
end

if ~isempty(img_idx) && ( min(img_idx) < 1 || max(img_idx) > hdr.dime.dim(5) )
  max_range = hdr.dime.dim(5);

  if max_range == 1
    error('"img_idx" should be 1.');
  else
    range = ['1 ' num2str(max_range)];
    error(['"img_idx" should be an integer within the range of [' range '].']);
  end
end

[img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,old_RGB);

return         % load_nii_img


%---------------------------------------------------------------------
function [img,hdr] = read_image(hdr, filetype,fileprefix,machine,img_idx,old_RGB)

switch filetype
  case {0, 1}
    fn = [fileprefix '.img'];
  case 2
    fn = [fileprefix '.nii'];
end

fid = fopen(fn,'r',machine);

if fid < 0,
  error('Cannot open file %s.',fn);
end

%  Set bitpix according to datatype
%
%  /*Acceptable values for datatype are*/
%
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64
%   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8
%   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%
switch hdr.dime.datatype
  case   1,
    hdr.dime.bitpix = 1;  precision = 'ubit1';
  case   2,
    hdr.dime.bitpix = 8;  precision = 'uint8';
  case   4,
    hdr.dime.bitpix = 16; precision = 'int16';
  case   8,
    hdr.dime.bitpix = 32; precision = 'int32';
  case  16,
    hdr.dime.bitpix = 32; precision = 'float32';
  case  32,
    hdr.dime.bitpix = 64; precision = 'float32';
  case  64,
    hdr.dime.bitpix = 64; precision = 'float64';
  case 128,
    hdr.dime.bitpix = 24; precision = 'uint8';
  case 256
    hdr.dime.bitpix = 8;  precision = 'int8';
  case 511
    hdr.dime.bitpix = 96; precision = 'float32';
  case 512
    hdr.dime.bitpix = 16; precision = 'uint16';
  case 768
    hdr.dime.bitpix = 32; precision = 'uint32';
  case 1024
    hdr.dime.bitpix = 64; precision = 'int64';
  case 1280
    hdr.dime.bitpix = 64; precision = 'uint64';
  case 1792,
    hdr.dime.bitpix = 128; precision = 'float64';
  otherwise
    error('This datatype is not supported');
end

if hdr.dime.dim(5) < 1
  hdr.dime.dim(5) = 1;
end

%  move pointer to the start of image block
%
switch filetype
  case {0, 1}
    fseek(fid, 0, 'bof');
  case 2
    fseek(fid, hdr.dime.vox_offset, 'bof');
end

%  Load whole image block for old Analyze format, or binary image,
%  or img_idx is empty; otherwise, load images that are specified
%  in img_idx
%
%  For binary image, we have to read all because pos can not be
%  seeked in bit and can not be calculated the way below.
%
if filetype == 0 || hdr.dime.datatype == 1 || isempty(img_idx)

  %  For each frame, precision of value will be read
  %  in img_siz times, where img_siz is only the
  %  dimension size of an image, not the byte storage
  %  size of an image.
  %
  img_siz = prod(hdr.dime.dim(2:5));

  %  For complex float32 or complex float64, voxel values
  %  include [real, imag]
  %
  if hdr.dime.datatype == 32 || hdr.dime.datatype == 1792
    img_siz = img_siz * 2;
  end

  %MPH: For RGB24, voxel values include 3 separate color planes
  %
  if hdr.dime.datatype == 128 || hdr.dime.datatype == 511
    img_siz = img_siz * 3;
  end

  img = fread(fid, img_siz, sprintf('*%s',precision));

else
  img = [];

  for i=1:length(img_idx)

    %  For each frame, precision of value will be read
    %  in img_siz times, where img_siz is only the
    %  dimension size of an image, not the byte storage
    %  size of an image.
    %
    img_siz = prod(hdr.dime.dim(2:4));

    %  Position is seeked in bytes. To convert dimension size
    %  to byte storage size, hdr.dime.bitpix/8 will be
    %  applied.
    %
    %  (MPH: This offset must be calculated BEFORE altering img_siz
    %  for the datatypes for which 'bitpix' and 'precision' do not
    %  "match").
    %
    pos = (img_idx(i) - 1) * img_siz * hdr.dime.bitpix/8;

    %  For complex float32 or complex float64, voxel values
    %  include [real, imag]
    %
    if hdr.dime.datatype == 32 || hdr.dime.datatype == 1792
      img_siz = img_siz * 2;
    end

    %MPH: For RGB24, voxel values include 3 separate color planes
    %
    if hdr.dime.datatype == 128 || hdr.dime.datatype == 511
      img_siz = img_siz * 3;
    end

    if filetype == 2
      fseek(fid, pos + hdr.dime.vox_offset, 'bof');
    else
      fseek(fid, pos, 'bof');
    end

    %  For each frame, fread will read precision of value
    %  in img_siz times
    %
    img = [img fread(fid, img_siz, sprintf('*%s',precision))]; %#ok
  end
end

%  For complex float32 or complex float64, voxel values
%  include [real, imag]
%
if hdr.dime.datatype == 32 || hdr.dime.datatype == 1792
  img = reshape(img, [2, length(img)/2]);
  img = complex(img(1,:)', img(2,:)');
end

fclose(fid);

%  Update the global min and max values
%
hdr.dime.glmax = max(double(img(:)));
hdr.dime.glmin = min(double(img(:)));

if isempty(img_idx)
  img_idx = 1:hdr.dime.dim(5);
end

if old_RGB && hdr.dime.datatype == 128 && hdr.dime.bitpix == 24
  img = squeeze(reshape(img, [hdr.dime.dim(2:3) 3 hdr.dime.dim(4) length(img_idx)]));
  img = permute(img, [1 2 4 3 5]);
elseif hdr.dime.datatype == 128 && hdr.dime.bitpix == 24
  img = squeeze(reshape(img, [3 hdr.dime.dim(2:4) length(img_idx)]));
  img = permute(img, [2 3 4 1 5]);
elseif hdr.dime.datatype == 511 && hdr.dime.bitpix == 96
  img = double(img);
  img = (img - min(img))/(max(img) - min(img));
  img = squeeze(reshape(img, [3 hdr.dime.dim(2:4) length(img_idx)]));
  img = permute(img, [2 3 4 1 5]);
else
  img = squeeze(reshape(img, [hdr.dime.dim(2:4) length(img_idx)]));
end

if ~isempty(img_idx)
  hdr.dime.dim(5) = length(img_idx);
end

return           % read_image


%---------------------------------------------------------------------
%  Load NIFTI dataset header. Support both *.nii and *.hdr/*.img file
%  extension. If file extension is not provided, *.hdr/*.img will be
%  used as default.
%
%  Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)
%
%  filename - NIFTI file name.
%
%  Returned values:
%
%  hdr - struct with NIFTI header fields.
%
%  filetype - 0 for Analyze format (*.hdr/*.img);
%     1 for NIFTI format in 2 files (*.hdr/*.img);
%     2 for NIFTI format in 1 file (*.nii).
%
%  fileprefix - NIFTI file name without extension.
%
%  machine    - a string, see below for details. The default here is 'ieee-le'.
%
%    'native'      or 'n' - local machine format - the default
%    'ieee-le'     or 'l' - IEEE floating point with little-endian
%                           byte ordering
%    'ieee-be'     or 'b' - IEEE floating point with big-endian
%                           byte ordering
%    'vaxd'        or 'd' - VAX D floating point and VAX ordering
%    'vaxg'        or 'g' - VAX G floating point and VAX ordering
%    'cray'        or 'c' - Cray floating point with big-endian
%                           byte ordering
%    'ieee-le.l64' or 'a' - IEEE floating point with little-endian
%                           byte ordering and 64 bit long data type
%    'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                           ordering and 64 bit long data type.
%
%  Number of scanned images in the file can be obtained by:
%  num_scan = hdr.dime.dim(5)
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%
%
% 2008/02/20 fixed error reading '*.nii.gz_restore.nii', replaced strrep

function [hdr, filetype, fileprefix, machine] = load_nii_hdr(fileprefix)

if ~exist('fileprefix','var'),
  error('Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)');
end

if ~exist('machine','var'), machine = 'ieee-le'; end

new_ext = 0;
[p,fname,ext]=fileparts(fileprefix);
if strcmpi('.nii',ext);
  new_ext = 1;
  fileprefix = fullfile(p,fname);
elseif strcmpi('.hdr', ext) || strcmpi('.img', ext)
  fileprefix = fullfile(p,fname);
end

%    if strfind('.nii',fileprefix)
%       new_ext = 1;
%       fileprefix = strrep(fileprefix,'.nii','');
%    end
%
%    if strfind('.hdr',fileprefix)
%       fileprefix = strrep(fileprefix,'.hdr','');
%    end
%
%    if strfind('.img',fileprefix)
%       fileprefix = strrep(fileprefix,'.img','');
%    end

if new_ext
  fn = sprintf('%s.nii',fileprefix);

  if ~exist(fn,'file')
    error('Cannot find file "%s.nii".', fileprefix);
  end
else
  fn = sprintf('%s.hdr',fileprefix);

  if ~exist(fn,'file')
    error('Cannot find file "%s.hdr".', fileprefix);
  end
end

fid = fopen(fn,'r',machine);

if fid < 0,
  error('Cannot open file %s.',fn);
else
  hdr = read_header(fid);
  fclose(fid);
end

if hdr.hk.sizeof_hdr ~= 348
  % first try reading the opposite endian to 'machine'
  switch machine,
    case 'ieee-le', machine = 'ieee-be';
    case 'ieee-be', machine = 'ieee-le';
  end

  fid = fopen(fn,'r',machine);

  if fid < 0,
    error('Cannot open file %s.',fn);
  else
    hdr = read_header(fid);
    fclose(fid);
  end
end

if hdr.hk.sizeof_hdr ~= 348
  % Now throw an error
  error('File "%s" is corrupted.',fn);
end

if strcmp(hdr.hist.magic, 'n+1')
  filetype = 2;
elseif strcmp(hdr.hist.magic, 'ni1')
  filetype = 1;
else
  filetype = 0;
end

return         % load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header(fid)

%  Original header structures
%  struct dsr
%       {
%       struct header_key1 hk;            /*   0 +  40       */
%       struct image_dimension1 dime;     /*  40 + 108       */
%       struct data_history1 hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

dsr.hk   = header_key1(fid);
dsr.dime = image_dimension1(fid);
dsr.hist = data_history1(fid);

%  For Analyze data format
%
if ~strcmp(dsr.hist.magic, 'n+1') && ~strcmp(dsr.hist.magic, 'ni1')
  dsr.hist.qform_code = 0;
  dsr.hist.sform_code = 0;
end

return          % read_header


%---------------------------------------------------------------------
function [ hk ] = header_key1(fid)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key1                    /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */
%
% int sizeof_header   Should be 348.
% char regular        Must be 'r' to indicate that all images and
%                     volumes are the same size.

hk.sizeof_hdr    = fread(fid, 1,'int32')';  % should be 348!
hk.data_type     = deblank(fread(fid,10,'*char')');
hk.db_name       = deblank(fread(fid,18,'*char')');
hk.extents       = fread(fid, 1,'int32')';
hk.session_error = fread(fid, 1,'int16')';
hk.regular       = fread(fid, 1,'*char')';
hk.dim_info      = fread(fid, 1,'char')';

return          % header_key1


%---------------------------------------------------------------------
function [ dime ] = image_dimension1(fid)

%  Original header structures
%  struct image_dimension1
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       /*
%           dim[0]      Number of dimensions in database; usually 4.
%           dim[1]      Image X dimension;  number of *pixels* in an image row.
%           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
%           dim[3]      Volume Z dimension; number of *slices* in a volume.
%           dim[4]      Time points; number of volumes in database
%       */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
% /*
%   pixdim[] specifies the voxel dimensions:
%   pixdim[1] - voxel width, mm
%   pixdim[2] - voxel height, mm
%   pixdim[3] - slice thickness, mm
%   pixdim[4] - volume timing, in msec
%         ..etc
% */
%       float vox_offset;                /* 68 + 4          */
%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
%       float scl_inter;   % float funused1;      /* 76 + 4 */
%       short slice_end;   % float funused2;      /* 80 + 2 */
%       char slice_code;   % float funused2;      /* 82 + 1 */
%       char xyzt_units;   % float funused2;      /* 83 + 1 */
%       float cal_max;                   /* 84 + 4          */
%       float cal_min;                   /* 88 + 4          */
%       float slice_duration;   % int compressed; /* 92 + 4 */
%       float toffset;   % int verified;          /* 96 + 4 */
%       int glmax;                       /* 100 + 4         */
%       int glmin;                       /* 104 + 4         */
%       };                               /* total=108 bytes */

dime.dim        = fread(fid,8,'int16')';
dime.intent_p1  = fread(fid,1,'float32')';
dime.intent_p2  = fread(fid,1,'float32')';
dime.intent_p3  = fread(fid,1,'float32')';
dime.intent_code = fread(fid,1,'int16')';
dime.datatype   = fread(fid,1,'int16')';
dime.bitpix     = fread(fid,1,'int16')';
dime.slice_start = fread(fid,1,'int16')';
dime.pixdim     = fread(fid,8,'float32')';
dime.vox_offset = fread(fid,1,'float32')';
dime.scl_slope  = fread(fid,1,'float32')';
dime.scl_inter  = fread(fid,1,'float32')';
dime.slice_end  = fread(fid,1,'int16')';
dime.slice_code = fread(fid,1,'char')';
dime.xyzt_units = fread(fid,1,'char')';
dime.cal_max    = fread(fid,1,'float32')';
dime.cal_min    = fread(fid,1,'float32')';
dime.slice_duration = fread(fid,1,'float32')';
dime.toffset    = fread(fid,1,'float32')';
dime.glmax      = fread(fid,1,'int32')';
dime.glmin      = fread(fid,1,'int32')';

return          % image_dimension1


%---------------------------------------------------------------------
function [ hist ] = data_history1(fid)

%  Original header structures
%  struct data_history1
%       {                                /* off + size      */
%       char descrip[80];                /* 0 + 80          */
%       char aux_file[24];               /* 80 + 24         */
%       short int qform_code;            /* 104 + 2         */
%       short int sform_code;            /* 106 + 2         */
%       float quatern_b;                 /* 108 + 4         */
%       float quatern_c;                 /* 112 + 4         */
%       float quatern_d;                 /* 116 + 4         */
%       float qoffset_x;                 /* 120 + 4         */
%       float qoffset_y;                 /* 124 + 4         */
%       float qoffset_z;                 /* 128 + 4         */
%       float srow_x[4];                 /* 132 + 16        */
%       float srow_y[4];                 /* 148 + 16        */
%       float srow_z[4];                 /* 164 + 16        */
%       char intent_name[16];            /* 180 + 16        */
%       char magic[4];   % int smin;     /* 196 + 4         */
%       };                               /* total=200 bytes */

hist.descrip     = deblank(fread(fid,80,'*char')');
hist.aux_file    = deblank(fread(fid,24,'*char')');
hist.qform_code  = fread(fid,1,'int16')';
hist.sform_code  = fread(fid,1,'int16')';
hist.quatern_b   = fread(fid,1,'float32')';
hist.quatern_c   = fread(fid,1,'float32')';
hist.quatern_d   = fread(fid,1,'float32')';
hist.qoffset_x   = fread(fid,1,'float32')';
hist.qoffset_y   = fread(fid,1,'float32')';
hist.qoffset_z   = fread(fid,1,'float32')';
hist.srow_x      = fread(fid,4,'float32')';
hist.srow_y      = fread(fid,4,'float32')';
hist.srow_z      = fread(fid,4,'float32')';
hist.intent_name = deblank(fread(fid,16,'*char')');
hist.magic       = deblank(fread(fid,4,'*char')');

fseek(fid,253,'bof');
hist.originator  = fread(fid, 5,'int16')';

return          % data_history1


%---------------------------------------------------------------------
function nii = load_niigz(fn, img_idx, old_RGB)
% nii = load_niigz(fn)
% load gzipped nifti file *.nii.gz
% see load_nii
% 2008/01/08 yamashiro, 2010/1027, mod for mac

if ~exist('img_idx','var'), img_idx = []; end
if ~exist('old_RGB','var'), old_RGB = 0; end

[p,f,e]=fileparts(fn);
if strcmpi(e,'.gz')
  wd = pwd;
  cd(p);

  %niifn = fullfile(p,f);
  niifn = ['tmp_' f];
  % -cd to avoid deteting original

  if isunix
    [s,w]=unix(sprintf('gzip -cd %s > %s',f,niifn));
  elseif ispc
    try
      [s,w]=dos(sprintf('gzip -cd %s > %s',f,niifn));
    catch
      [s,w]=dos(sprintf('"%s%sbin%sgzip.exe " -cd "%s" > "%s"',fileparts(mfilename('fullpath')),filesep(),filesep(),f,niifn));
    end
  end
  if s~=0
    error(w);
  end
  nii = load_nii(niifn, img_idx, old_RGB);
  nii.fileprefix = nii.fileprefix(5:end); % remove tmp_
  eval(sprintf('!rm %s',niifn));
  cd(wd);
else % not .gz
  if strcmpi(e,'.nii')
    nii = load_nii(fn, img_idx, old_RGB);
  else
    error(['cannot read ' fn]);
  end
end

%     %niifn = fullfile(p,f);
%     niifn = fullfile(p,['tmp_' f])
%     % -cd to avoid deteting original
%     %eval(sprintf('!gzip -cd %s > %s', fn, niifn));
%
%     if isunix
%         [s,w]=unix(sprintf('gzip -cd %s > %s',f,niifn));
%     elseif ispc
%         try
%             [s,w]=dos(sprintf('gzip -cd %s > %s',f,niifn));
%         catch
%             [s,w]=dos(sprintf('"%s%sbin%sgzip.exe " -cd "%s" > "%s"',fileparts(mfilename('fullpath')),filesep(),filesep()),f,niifn));
%         end
%     end
%     if s~=0
%         error(w);
%     end
%     nii = load_nii(niifn, img_idx, old_RGB);
%     nii.fileprefix = nii.fileprefix(5:end); % remove tmp_
%     eval(sprintf('!rm %s',niifn));


%---------------------------------------------------------------------
%  Perform a subset of NIFTI sform/qform transform. Transforms like
%  (Translation, Flipping, and a few Rotation (N*90 degree) are
%  supported. Other transforms (any degree rotation, shears, etc.)
%  are not supported, because in those transforms, each voxel has
%  to be repositioned, interpolated, and whole image(s) will have
%  to be reconstructed. If an input data (nii) can not be handled,
%  The program will exit with an error message "Transform of this
%  NIFTI data is not supported by the program". After the transform,
%  nii will be in RAS orientation, i.e. X axis from Left to Right,
%  Y axis from Posterior to Anterior, and Z axis from Inferior to
%  Superior. The RAS orientation system sometimes is also referred
%  as right-hand coordinate system, or Neurologist preferred system.
%
%  NOTE: This function should be called immediately after load_nii,
%        if loaded data will be sent to view_nii, plsgui, etc.
%
%  Usage: [ nii ] = xform_nii(nii)
%
%  nii  - NIFTI structure (returned from load_nii)
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  7-nov-2005: Images whose cardinal planes are slightly off Cartesian
% coordinates will also be loaded. However, this approximation
% is based on the following assumption: In affine matrix, any
% value (absolute) below a tenth of the third largest value
% (absolute) can be ignored and replaced with 0. In this case,
% fields 'old_affine' & 'new_affine' will be added into hdr.hist.
% If you can not accept the above assumption, simply discard any
% nii structure with fields 'old_affine' & 'new_affine'.
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%
function nii = xform_nii(nii)

%  save a copy of the header as it was loaded.  This is the
%  header before any sform, qform manipulation is done.
%
nii.original.hdr = nii.hdr;

%  if scl_slope field is nonzero, then each voxel value in the
%  dataset should be scaled as: y = scl_slope * x + scl_inter
%  I bring it here because hdr will be modified by change_hdr.
%
if nii.hdr.dime.scl_slope ~= 0 && ismember(nii.hdr.dime.datatype, [2,4,8,16,64,256,512,768])

  nii.img = ...
      nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;
  nii.hdr.dime.datatype = 64;
  nii.hdr.dime.bitpix = 64;
  nii.hdr.dime.glmax = max(nii.img(:));
  nii.hdr.dime.glmin = min(nii.img(:));

  %  set scale to non-use, because it is applied in xform_nii
  %
  nii.hdr.dime.scl_slope = 0;
end

%  However, the scaling is to be ignored if datatype is DT_RGB24.

%  If datatype is a complex type, then the scaling is to be applied
%  to both the real and imaginary parts.
%
if nii.hdr.dime.scl_slope ~= 0 && ismember(nii.hdr.dime.datatype, [32,1792])

  nii.img = ...
      nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;
  if nii.hdr.dime.datatype == 32
    nii.img = single(nii.img);
  end

  %  set scale to non-use, because it is applied in xform_nii
  %
  nii.hdr.dime.scl_slope = 0;
end

%  There is no need for this program to transform Analyze data
%
if nii.filetype == 0
  nii.hdr.hist.rot_orient = [];
  nii.hdr.hist.flip_orient = [];

  return;     % no sform/qform for Analyze format
end

hdr = nii.hdr;

[hdr, orient] = change_hdr(hdr);

%  flip and/or rotate image data
%
if ~isequal(orient, [1 2 3])

  old_dim = hdr.dime.dim(2:4);

  %  More than 1 time frame
  %
  if ndims(nii.img) > 3
    pattern = 1:prod(old_dim);
  else
    pattern = [];
  end

  if ~isempty(pattern)
    pattern = reshape(pattern, old_dim);
  end

  %  calculate for rotation after flip
  %
  rot_orient = mod(orient + 2, 3) + 1;

  %  do flip:
  %
  flip_orient = orient - rot_orient;

  for i = 1:3
    if flip_orient(i)
      if ~isempty(pattern)
        pattern = flipdim(pattern, i); %#ok
      else
        nii.img = flipdim(nii.img, i); %#ok
      end
    end
  end

  %  get index of orient (rotate inversely)
  %
  [tmp, rot_orient] = sort(rot_orient); %#ok

  new_dim = old_dim;
  new_dim = new_dim(rot_orient);
  hdr.dime.dim(2:4) = new_dim;

  new_pixdim = hdr.dime.pixdim(2:4);
  new_pixdim = new_pixdim(rot_orient);
  hdr.dime.pixdim(2:4) = new_pixdim;

  %  re-calculate originator
  %
  tmp = hdr.hist.originator(1:3);
  tmp = tmp(rot_orient);
  flip_orient = flip_orient(rot_orient);

  for i = 1:3
    if flip_orient(i) && ~isequal(tmp(i), 0)
      tmp(i) = new_dim(i) - tmp(i) + 1;
    end
  end

  hdr.hist.originator(1:3) = tmp;
  hdr.hist.rot_orient = rot_orient;
  hdr.hist.flip_orient = flip_orient;

  %  do rotation:
  %
  if ~isempty(pattern)
    pattern = permute(pattern, rot_orient);
    pattern = pattern(:);

    nii.img = reshape(nii.img, [prod(new_dim) hdr.dime.dim(5)]);
    nii.img = nii.img(pattern, :);
    nii.img = reshape(nii.img, [new_dim       hdr.dime.dim(5)]);
  else
    nii.img = permute(nii.img, rot_orient);
  end
else
  hdr.hist.rot_orient = [];
  hdr.hist.flip_orient = [];
end

nii.hdr = hdr;

return;          % xform_nii


%-----------------------------------------------------------------------
function [hdr, orient] = change_hdr(hdr)

orient = [1 2 3];
affine_transform = 1;

%  NIFTI can have both sform and qform transform. This program
%  will check sform_code prior to qform_code
%
if hdr.hist.sform_code > 0
  R = [hdr.hist.srow_x(1:3)
       hdr.hist.srow_y(1:3)
       hdr.hist.srow_z(1:3)];

  T = [hdr.hist.srow_x(4)
       hdr.hist.srow_y(4)
       hdr.hist.srow_z(4)];
elseif hdr.hist.qform_code > 0
  b = hdr.hist.quatern_b;
  c = hdr.hist.quatern_c;
  d = hdr.hist.quatern_d;

  if 1.0-(b*b+c*c+d*d) < 0
    if abs(1.0-(b*b+c*c+d*d)) < 1e-5
      a = 0;
    else
      error('Incorrect quaternion values in this NIFTI data.');
    end
  else
    a = sqrt(1.0-(b*b+c*c+d*d));
  end

  qfac = hdr.dime.pixdim(1);
  i = hdr.dime.pixdim(2);
  j = hdr.dime.pixdim(3);
  k = qfac * hdr.dime.pixdim(4);

  R = [a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
       2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
       2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];

  R = R * diag([i j k]);

  T = [hdr.hist.qoffset_x
       hdr.hist.qoffset_y
       hdr.hist.qoffset_z];
else
  affine_transform = 0; % no sform or qform transform
end

if affine_transform == 1
  if det(R) == 0 || ~isequal(R(find(R)), sum(R)') %#ok
    hdr.hist.old_affine = R;
    R_sort = sort(abs(R(:)));
    R( find( abs(R) < min(R_sort(end-2:end))/10 ) ) = 0; %#ok
    hdr.hist.new_affine = R;

    if det(R) == 0 || ~isequal(R(find(R)), sum(R)') %#ok
      error('Transform of this NIFTI data is not supported by the program.');
    end
  end

  voxel_size = abs(sum(R,1));
  inv_R = inv(R);
  originator = round(abs(inv_R*T)+1); %#ok
  orient = get_orient(inv_R);

  %  modify pixdim and originator
  %
  hdr.dime.pixdim(2:4) = voxel_size;
  hdr.hist.originator(1:3) = originator;

  %  set sform or qform to non-use, because they have been
  %  applied in xform_nii
  %
  hdr.hist.qform_code = 0;
  hdr.hist.sform_code = 0;
end

%  apply space_unit to pixdim if not 1 (mm)
%
space_unit = get_units(hdr);

if space_unit ~= 1
  hdr.dime.pixdim(2:4) = hdr.dime.pixdim(2:4) * space_unit;

  %  set space_unit of xyzt_units to millimeter, because
  %  voxel_size has been re-scaled
  %
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,1,0));
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,2,1));
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,3,0));
end

return;          % change_hdr


%-----------------------------------------------------------------------
function orient = get_orient(R)

orient = [];

for i = 1:3
  switch find(R(i,:)) * sign(sum(R(i,:)))
    case 1
      orient = [orient 1];   %#ok % Left to Right
    case 2
      orient = [orient 2];   %#ok % Posterior to Anterior
    case 3
      orient = [orient 3];   %#ok % Inferior to Superior
    case -1
      orient = [orient 4];   %#ok % Right to Left
    case -2
      orient = [orient 5];   %#ok % Anterior to Posterior
    case -3
      orient = [orient 6];   %#ok % Superior to Inferior
  end
end

return;          % get_orient


%-----------------------------------------------------------------------
function [space_unit, time_unit] = get_units(hdr)

switch bitand(hdr.dime.xyzt_units, 7)  % mask with 0x07
  case 1
    space_unit = 1e+3;    % meter, m
  case 3
    space_unit = 1e-3;    % micrometer, um
  otherwise
    space_unit = 1;     % millimeter, mm
end

switch bitand(hdr.dime.xyzt_units, 56) % mask with 0x38
  case 16
    time_unit = 1e-3;     % millisecond, ms
  case 24
    time_unit = 1e-6;     % microsecond, us
  otherwise
    time_unit = 1;      % second, s
end

return;          % get_units


%-----------------------------------------------------------------------
%  Save NIFTI dataset. Support both *.nii and *.hdr/*.img file extension.
%  If file extension is not provided, *.hdr/*.img will be used as default.
%
%  Usage: save_nii(nii, filename, [old_RGB])
%
%  nii.hdr - struct with NIFTI header fields.
%  nii.img - 3D (or 4D) matrix of NIFTI data.
%  filename - NIFTI file name.
%
%  old_RGB    - an optional boolean variable to handle special RGB data
%       sequence [R1 R2 ... G1 G2 ... B1 B2 ...] that is used only by
%       AnalyzeDirect (Analyze Software). Since both NIfTI and Analyze
%       file format use RGB triple [R1 G1 B1 R2 G2 B2 ...] sequentially
%       for each voxel, this variable is set to FALSE by default. If you
%       would like the saved image only to be opened by AnalyzeDirect
%       Software, set old_RGB to TRUE (or 1).
%
%  Tip: to change the data type, set nii.hdr.dime.datatype,
% and nii.hdr.dime.bitpix to:
%
%     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN
%     1 Binary                         (ubit1, bitpix=1) % DT_BINARY
%     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8
%     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32
%    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32
%    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
%    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64
%   128 Red-Green-Blue            (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24
%   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8
%   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16
%   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32
%  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
%  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64
%  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128
%  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%  - "old_RGB" related codes in "save_nii.m" are added by Mike Harms (2006.06.28)
%
function save_nii(nii, fileprefix, old_RGB)

if ~exist('nii','var') || isempty(nii) || ~isfield(nii,'hdr') || ...
      ~isfield(nii,'img') || ~exist('fileprefix','var') || isempty(fileprefix)

  error('Usage: save_nii(nii, filename, [old_RGB])');
end

if ~exist('old_RGB','var'), old_RGB = 0; end

filetype = 1;

if strfind('.nii',fileprefix) %#ok
  filetype = 2;
  fileprefix = strrep(fileprefix,'.nii','');
end

if strfind('.hdr',fileprefix) %#ok
  fileprefix = strrep(fileprefix,'.hdr','');
end

if strfind('.img',fileprefix) %#ok
  fileprefix = strrep(fileprefix,'.img','');
end

write_nii(nii, filetype, fileprefix, old_RGB);

return          % save_nii


%-----------------------------------------------------------------------------------
function write_nii(nii, filetype, fileprefix, old_RGB)

hdr = nii.hdr;

switch double(hdr.dime.datatype),
  case   1,
    hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
  case   2,
    hdr.dime.bitpix = int16(8 ); precision = 'uint8';
  case   4,
    hdr.dime.bitpix = int16(16); precision = 'int16';
  case   8,
    hdr.dime.bitpix = int16(32); precision = 'int32';
  case  16,
    hdr.dime.bitpix = int16(32); precision = 'float32';
  case  32,
    hdr.dime.bitpix = int16(64); precision = 'float32';
  case  64,
    hdr.dime.bitpix = int16(64); precision = 'float64';
  case 128,
    hdr.dime.bitpix = int16(24); precision = 'uint8';
  case 256
    hdr.dime.bitpix = int16(8 ); precision = 'int8';
  case 512
    hdr.dime.bitpix = int16(16); precision = 'uint16';
  case 768
    hdr.dime.bitpix = int16(32); precision = 'uint32';
  case 1024
    hdr.dime.bitpix = int16(64); precision = 'int64';
  case 1280
    hdr.dime.bitpix = int16(64); precision = 'uint64';
  case 1792,
    hdr.dime.bitpix = int16(128); precision = 'float64';
  otherwise
    error('This datatype is not supported');
end

hdr.dime.glmax = round(double(max(nii.img(:))));
hdr.dime.glmin = round(double(min(nii.img(:))));

if filetype == 2
  fid = fopen(sprintf('%s.nii',fileprefix),'w');

  if fid < 0,
    msg = sprintf('Cannot open file %s.nii.',fileprefix);
    error(msg); %#ok
  end

  hdr.dime.vox_offset = 352;
  hdr.hist.magic = 'n+1';
  save_nii_hdr(hdr, fid);
else
  fid = fopen(sprintf('%s.hdr',fileprefix),'w');

  if fid < 0,
    msg = sprintf('Cannot open file %s.hdr.',fileprefix);
    error(msg); %#ok
  end

  hdr.dime.vox_offset = 0;
  hdr.hist.magic = 'ni1';
  save_nii_hdr(hdr, fid);

  fclose(fid);
  fid = fopen(sprintf('%s.img',fileprefix),'w');
end

ScanDim = double(hdr.dime.dim(5));     %#ok % t
SliceDim = double(hdr.dime.dim(4));    %#ok % z
RowDim   = double(hdr.dime.dim(3));    %#ok % y
PixelDim = double(hdr.dime.dim(2));         % x
SliceSz  = double(hdr.dime.pixdim(4)); %#ok
RowSz    = double(hdr.dime.pixdim(3)); %#ok
PixelSz  = double(hdr.dime.pixdim(2)); %#ok

x = 1:PixelDim; %#ok

if filetype == 2
  skip_bytes = double(hdr.dime.vox_offset) - 348;
else
  skip_bytes = 0;
end

if double(hdr.dime.datatype) == 128

  %  RGB planes are expected to be in the 4th dimension of nii.img
  %
  if(size(nii.img,4)~=3)
    error('The NII structure does not appear to have 3 RGB color planes in the 4th dimension');
  end

  if old_RGB
    nii.img = permute(nii.img, [1 2 4 3 5]);
  else
    nii.img = permute(nii.img, [4 1 2 3 5]);
  end
end

%  For complex float32 or complex float64, voxel values
%  include [real, imag]
%
if hdr.dime.datatype == 32 || hdr.dime.datatype == 1792
  real_img = real(nii.img(:))';
  nii.img = imag(nii.img(:))';
  nii.img = [real_img; nii.img];
end

if skip_bytes
  fwrite(fid, ones(1,skip_bytes), 'uint8');
end

fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
fclose(fid);

return;         % write_nii


%-----------------------------------------------------------------------------------
function save_niigz(nii, fileprefix, old_RGB)
% save_niigz(nii, fn)
% save nifti file in *.nii.gz format
% see save_nii
% 2007/01/08 yamashiro

if ~exist('old_RGB','var'), old_RGB = 0; end
if strfind('.nii.gz',fileprefix) %#ok % .nii.gz
  niifn = strrep(fileprefix,'.nii.gz', '.nii');
elseif strfind('.nii',fileprefix) %#ok % .nii
  niifn = fileprefix;
else % no ext
  niifn = [fileprefix, '.nii'];
end
gzfn = [niifn, '.gz']; %#ok
save_nii(nii,niifn,old_RGB);

if isunix
  eval(sprintf('!gzip %s',niifn));
elseif ispc
  try
    eval(sprintf('!gzip %s',niifn));
  catch
    eval(sprintf('!"%s%sbin%sgzip.exe " "%s"',fileparts(mfilename('fullpath')),filesep(),filesep(),niifn));
  end
end

return

if exist(gzfn,'file'); %#ok
  fprintf('%s will be over written\n',gzfn);
  delete(gzfn);
end

if isunix
  [s,w]=unix(['gzip ',niifn]);
elseif ispc
  try
    [s,w]=dos(['gzip ',niifn]);
  catch
    [s,w]=dos([sprintf('"%s%sbin%sgzip.exe " ',fileparts(mfilename('fullpath')),filesep(),filesep()),'"',niifn,'"']);
  end
end

if s==0
  fprintf('wrote %s\n', gzfn);
else
  error(w);
end
return


%-----------------------------------------------------------------------------------
%  Save NIFTI dataset header. Support both *.nii and *.hdr/*.img file
%  extension.
%
%  Usage: save_nii_hdr(hdr, fid)
%
%  hdr - struct with NIFTI header fields.
%
%  fileprefix - NIFTI file name without extension.
%
%  Part of this file is copied and modified under GNU license from
%  MRI_TOOLBOX developed by CNSP in Flinders University, Australia
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%
function save_nii_hdr(hdr, fid)

if ~exist('hdr','var') || ~exist('fid','var')
  error('Usage: save_nii_hdr(hdr, fid)');
end

if ~isequal(hdr.hk.sizeof_hdr,348),
  error('hdr.hk.sizeof_hdr must be 348.');
end

if hdr.hist.qform_code == 0 && hdr.hist.sform_code == 0
  hdr.hist.sform_code = 1;
  hdr.hist.srow_x(1) = hdr.dime.pixdim(2);
  hdr.hist.srow_x(2) = 0;
  hdr.hist.srow_x(3) = 0;
  hdr.hist.srow_y(1) = 0;
  hdr.hist.srow_y(2) = hdr.dime.pixdim(3);
  hdr.hist.srow_y(3) = 0;
  hdr.hist.srow_z(1) = 0;
  hdr.hist.srow_z(2) = 0;
  hdr.hist.srow_z(3) = hdr.dime.pixdim(4);
  hdr.hist.srow_x(4) = (1-hdr.hist.originator(1))*hdr.dime.pixdim(2);
  hdr.hist.srow_y(4) = (1-hdr.hist.originator(2))*hdr.dime.pixdim(3);
  hdr.hist.srow_z(4) = (1-hdr.hist.originator(3))*hdr.dime.pixdim(4);
end

write_header(hdr, fid);

return;         % save_nii_hdr


%---------------------------------------------------------------------
function write_header(hdr, fid)

%  Original header structures
%  struct dsr       /* dsr = hdr */
%       {
%       struct header_key2 hk;            /*   0 +  40       */
%       struct image_dimension2 dime;     /*  40 + 108       */
%       struct data_history2 hist;        /* 148 + 200       */
%       };                               /* total= 348 bytes*/

header_key2(fid, hdr.hk);
image_dimension2(fid, hdr.dime);
data_history2(fid, hdr.hist);

%  check the file size is 348 bytes
%
fbytes = ftell(fid);

if ~isequal(fbytes,348),
  msg = sprintf('Header size is not 348 bytes.');
  warning(msg); %#ok
end

return;         % write_header


%---------------------------------------------------------------------
function header_key2(fid, hk)

fseek(fid,0,'bof');

%  Original header structures
%  struct header_key2                      /* header key      */
%       {                                /* off + size      */
%       int sizeof_hdr                   /*  0 +  4         */
%       char data_type[10];              /*  4 + 10         */
%       char db_name[18];                /* 14 + 18         */
%       int extents;                     /* 32 +  4         */
%       short int session_error;         /* 36 +  2         */
%       char regular;                    /* 38 +  1         */
%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
%       };                               /* total=40 bytes  */

fwrite(fid, hk.sizeof_hdr(1),    'int32');  % must be 348.

% data_type = sprintf('%-10s',hk.data_type);  % ensure it is 10 chars from left
% fwrite(fid, data_type(1:10), 'uchar');
pad = zeros(1, 10-length(hk.data_type));
hk.data_type = [hk.data_type  char(pad)];
fwrite(fid, hk.data_type(1:10), 'uchar');

% db_name   = sprintf('%-18s', hk.db_name); % ensure it is 18 chars from left
% fwrite(fid, db_name(1:18), 'uchar');
pad = zeros(1, 18-length(hk.db_name));
hk.db_name = [hk.db_name  char(pad)];
fwrite(fid, hk.db_name(1:18), 'uchar');

fwrite(fid, hk.extents(1),       'int32');
fwrite(fid, hk.session_error(1), 'int16');
fwrite(fid, hk.regular(1),       'uchar');  % might be uint8

% fwrite(fid, hk.hkey_un0(1),    'uchar');
% fwrite(fid, hk.hkey_un0(1),    'uint8');
fwrite(fid, hk.dim_info(1),      'uchar');

return;         % header_key2


%---------------------------------------------------------------------
function image_dimension2(fid, dime)

%  Original header structures
%  struct image_dimension2
%       {                                /* off + size      */
%       short int dim[8];                /* 0 + 16          */
%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
%       short int intent_code;   % short int unused1;   /* 28 + 2 */
%       short int datatype;              /* 30 + 2          */
%       short int bitpix;                /* 32 + 2          */
%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
%       float pixdim[8];                 /* 36 + 32         */
%     /*
%       pixdim[] specifies the voxel dimensions:
%       pixdim[1] - voxel width
%       pixdim[2] - voxel height
%       pixdim[3] - interslice distance
%       pixdim[4] - volume timing, in msec
%         ..etc
%     */
%       float vox_offset;                /* 68 + 4          */
%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
%       float scl_inter;   % float funused1;      /* 76 + 4 */
%       short slice_end;   % float funused2;      /* 80 + 2 */
%       char slice_code;   % float funused2;      /* 82 + 1 */
%       char xyzt_units;   % float funused2;      /* 83 + 1 */
%       float cal_max;                   /* 84 + 4          */
%       float cal_min;                   /* 88 + 4          */
%       float slice_duration;   % int compressed; /* 92 + 4 */
%       float toffset;   % int verified;          /* 96 + 4 */
%       int glmax;                       /* 100 + 4         */
%       int glmin;                       /* 104 + 4         */
%       };                               /* total=108 bytes */

fwrite(fid, dime.dim(1:8),        'int16');
fwrite(fid, dime.intent_p1(1),  'float32');
fwrite(fid, dime.intent_p2(1),  'float32');
fwrite(fid, dime.intent_p3(1),  'float32');
fwrite(fid, dime.intent_code(1),  'int16');
fwrite(fid, dime.datatype(1),     'int16');
fwrite(fid, dime.bitpix(1),       'int16');
fwrite(fid, dime.slice_start(1),  'int16');
fwrite(fid, dime.pixdim(1:8),   'float32');
fwrite(fid, dime.vox_offset(1), 'float32');
fwrite(fid, dime.scl_slope(1),  'float32');
fwrite(fid, dime.scl_inter(1),  'float32');
fwrite(fid, dime.slice_end(1),    'int16');
fwrite(fid, dime.slice_code(1),   'uchar');
fwrite(fid, dime.xyzt_units(1),   'uchar');
fwrite(fid, dime.cal_max(1),    'float32');
fwrite(fid, dime.cal_min(1),    'float32');
fwrite(fid, dime.slice_duration(1), 'float32');
fwrite(fid, dime.toffset(1),    'float32');
fwrite(fid, dime.glmax(1),        'int32');
fwrite(fid, dime.glmin(1),        'int32');

return;         % image_dimension2


%---------------------------------------------------------------------
function data_history2(fid, hist)

% Original header structures
%struct data_history2
%       {                                /* off + size      */
%       char descrip[80];                /* 0 + 80          */
%       char aux_file[24];               /* 80 + 24         */
%       short int qform_code;            /* 104 + 2         */
%       short int sform_code;            /* 106 + 2         */
%       float quatern_b;                 /* 108 + 4         */
%       float quatern_c;                 /* 112 + 4         */
%       float quatern_d;                 /* 116 + 4         */
%       float qoffset_x;                 /* 120 + 4         */
%       float qoffset_y;                 /* 124 + 4         */
%       float qoffset_z;                 /* 128 + 4         */
%       float srow_x[4];                 /* 132 + 16        */
%       float srow_y[4];                 /* 148 + 16        */
%       float srow_z[4];                 /* 164 + 16        */
%       char intent_name[16];            /* 180 + 16        */
%       char magic[4];   % int smin;     /* 196 + 4         */
%       };                               /* total=200 bytes */

% descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
% fwrite(fid, descrip(1:80),    'uchar');
pad = zeros(1, 80-length(hist.descrip));
hist.descrip = [hist.descrip  char(pad)];
fwrite(fid, hist.descrip(1:80), 'uchar');

% aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
% fwrite(fid, aux_file(1:24),   'uchar');
pad = zeros(1, 24-length(hist.aux_file));
hist.aux_file = [hist.aux_file  char(pad)];
fwrite(fid, hist.aux_file(1:24), 'uchar');

fwrite(fid, hist.qform_code,    'int16');
fwrite(fid, hist.sform_code,    'int16');
fwrite(fid, hist.quatern_b,   'float32');
fwrite(fid, hist.quatern_c,   'float32');
fwrite(fid, hist.quatern_d,   'float32');
fwrite(fid, hist.qoffset_x,   'float32');
fwrite(fid, hist.qoffset_y,   'float32');
fwrite(fid, hist.qoffset_z,   'float32');
fwrite(fid, hist.srow_x(1:4), 'float32');
fwrite(fid, hist.srow_y(1:4), 'float32');
fwrite(fid, hist.srow_z(1:4), 'float32');

% intent_name = sprintf('%-16s', hist.intent_name); % 16 chars from left
% fwrite(fid, intent_name(1:16),    'uchar');
pad = zeros(1, 16-length(hist.intent_name));
hist.intent_name = [hist.intent_name  char(pad)];
fwrite(fid, hist.intent_name(1:16), 'uchar');

% magic = sprintf('%-4s', hist.magic);    % 4 chars from left
% fwrite(fid, magic(1:4),           'uchar');
pad = zeros(1, 4-length(hist.magic));
hist.magic = [hist.magic  char(pad)];
fwrite(fid, hist.magic(1:4),        'uchar');

return;         % data_history2

%  Perform a subset of NIFTI sform/qform transform. Transforms like
%  (Translation, Flipping, and a few Rotation (N*90 degree) are
%  supported. Other transforms (any degree rotation, shears, etc.)
%  are not supported, because in those transforms, each voxel has
%  to be repositioned, interpolated, and whole image(s) will have
%  to be reconstructed. If an input data (nii) can not be handled,
%  The program will exit with an error message "Transform of this
%  NIFTI data is not supported by the program". After the transform,
%  nii will be in RAS orientation, i.e. X axis from Left to Right,
%  Y axis from Posterior to Anterior, and Z axis from Inferior to
%  Superior. The RAS orientation system sometimes is also referred
%  as right-hand coordinate system, or Neurologist preferred system.
%
%  NOTE: This function should be called immediately after load_nii,
%        if loaded data will be sent to view_nii, plsgui, etc.
%
%  Usage: [ nii ] = xform_nii(nii)
%
%  nii  - NIFTI structure (returned from load_nii)
%
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%
%  7-nov-2005: Images whose cardinal planes are slightly off Cartesian
% coordinates will also be loaded. However, this approximation
% is based on the following assumption: In affine matrix, any
% value (absolute) below a tenth of the third largest value
% (absolute) can be ignored and replaced with 0. In this case,
% fields 'old_affine' & 'new_affine' will be added into hdr.hist.
% If you can not accept the above assumption, simply discard any
% nii structure with fields 'old_affine' & 'new_affine'.
%
%  - Jimmy Shen (pls@rotman-baycrest.on.ca)
%
function nii = xform_nii(nii)

%  save a copy of the header as it was loaded.  This is the
%  header before any sform, qform manipulation is done.
%
nii.original.hdr = nii.hdr;

%  if scl_slope field is nonzero, then each voxel value in the
%  dataset should be scaled as: y = scl_slope * x + scl_inter
%  I bring it here because hdr will be modified by change_hdr.
%
if nii.hdr.dime.scl_slope ~= 0 && ...
      ismember(nii.hdr.dime.datatype, [2,4,8,16,64,256,512,768])

  nii.img = ...
      nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;
  nii.hdr.dime.datatype = 64;
  nii.hdr.dime.bitpix = 64;
  nii.hdr.dime.glmax = max(nii.img(:));
  nii.hdr.dime.glmin = min(nii.img(:));

  %  set scale to non-use, because it is applied in xform_nii
  %
  nii.hdr.dime.scl_slope = 0;
end

%  However, the scaling is to be ignored if datatype is DT_RGB24.

%  If datatype is a complex type, then the scaling is to be applied
%  to both the real and imaginary parts.
%
if nii.hdr.dime.scl_slope ~= 0 && ...
      ismember(nii.hdr.dime.datatype, [32,1792])

  nii.img = ...
      nii.hdr.dime.scl_slope * double(nii.img) + nii.hdr.dime.scl_inter;
  if nii.hdr.dime.datatype == 32
    nii.img = single(nii.img);
  end

  %  set scale to non-use, because it is applied in xform_nii
  %
  nii.hdr.dime.scl_slope = 0;
end

%  There is no need for this program to transform Analyze data
%
if nii.filetype == 0
  nii.hdr.hist.rot_orient = [];
  nii.hdr.hist.flip_orient = [];

  return;     % no sform/qform for Analyze format
end

hdr = nii.hdr;

[hdr,orient] = change_hdr(hdr);

%  flip and/or rotate image data
%
if ~isequal(orient, [1 2 3])

  old_dim = hdr.dime.dim(2:4);

  %  More than 1 time frame
  %
  if ndims(nii.img) > 3
    pattern = 1:prod(old_dim);
  else
    pattern = [];
  end

  if ~isempty(pattern)
    pattern = reshape(pattern, old_dim);
  end

  %  calculate for rotation after flip
  %
  rot_orient = mod(orient + 2, 3) + 1;

  %  do flip:
  %
  flip_orient = orient - rot_orient;

  for i = 1:3
    if flip_orient(i)
      if ~isempty(pattern)
        pattern = flip(pattern, i);
      else
        nii.img = flip(nii.img, i);
      end
    end
  end

  %  get index of orient (rotate inversely)
  %
  [tmp,rot_orient] = sort(rot_orient); %#ok

  new_dim = old_dim;
  new_dim = new_dim(rot_orient);
  hdr.dime.dim(2:4) = new_dim;

  new_pixdim = hdr.dime.pixdim(2:4);
  new_pixdim = new_pixdim(rot_orient);
  hdr.dime.pixdim(2:4) = new_pixdim;

  %  re-calculate originator
  %
  tmp = hdr.hist.originator(1:3);
  tmp = tmp(rot_orient);
  flip_orient = flip_orient(rot_orient);

  for i = 1:3
    if flip_orient(i) && ~isequal(tmp(i), 0)
      tmp(i) = new_dim(i) - tmp(i) + 1;
    end
  end

  hdr.hist.originator(1:3) = tmp;
  hdr.hist.rot_orient = rot_orient;
  hdr.hist.flip_orient = flip_orient;

  %  do rotation:
  %
  if ~isempty(pattern)
    pattern = permute(pattern, rot_orient);
    pattern = pattern(:);

    nii.img = reshape(nii.img, [prod(new_dim) hdr.dime.dim(5)]);
    nii.img = nii.img(pattern, :);
    nii.img = reshape(nii.img, [new_dim       hdr.dime.dim(5)]);
  else
    nii.img = permute(nii.img, rot_orient);
  end
else
  hdr.hist.rot_orient = [];
  hdr.hist.flip_orient = [];
end

nii.hdr = hdr;

return;         % xform_nii


%-----------------------------------------------------------------------
function [hdr, orient] = change_hdr(hdr)

orient = [1 2 3];
affine_transform = 1;

%  NIFTI can have both sform and qform transform. This program
%  will check sform_code prior to qform_code
%
if hdr.hist.sform_code > 0
  R = [hdr.hist.srow_x(1:3)
       hdr.hist.srow_y(1:3)
       hdr.hist.srow_z(1:3)];

  T = [hdr.hist.srow_x(4)
       hdr.hist.srow_y(4)
       hdr.hist.srow_z(4)];
elseif hdr.hist.qform_code > 0
  b = hdr.hist.quatern_b;
  c = hdr.hist.quatern_c;
  d = hdr.hist.quatern_d;

  if 1.0-(b*b+c*c+d*d) < 0
    if abs(1.0-(b*b+c*c+d*d)) < 1e-5
      a = 0;
    else
      error('Incorrect quaternion values in this NIFTI data.');
    end
  else
    a = sqrt(1.0-(b*b+c*c+d*d));
  end

  qfac = hdr.dime.pixdim(1);
  i = hdr.dime.pixdim(2);
  j = hdr.dime.pixdim(3);
  k = qfac * hdr.dime.pixdim(4);

  R = [a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
       2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
       2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];

  R = R * diag([i j k]);

  T = [hdr.hist.qoffset_x
       hdr.hist.qoffset_y
       hdr.hist.qoffset_z];
else
  affine_transform = 0; % no sform or qform transform
end

if affine_transform == 1
  if det(R) == 0 || ~isequal(R(find(R)), sum(R)') %#ok
    hdr.hist.old_affine = R;
    R_sort = sort(abs(R(:)));
    R( find( abs(R) < min(R_sort(end-2:end))/10 ) ) = 0; %#ok
    hdr.hist.new_affine = R;

    if det(R) == 0 || ~isequal(R(find(R)), sum(R)') %#ok
      error('Transform of this NIFTI data is not supported by the program.');
    end
  end

  voxel_size = abs(sum(R,1));
  inv_R = inv(R);
  originator = round(abs(inv_R*T)+1); %#ok
  orient = get_orient(inv_R);

  %  modify pixdim and originator
  %
  hdr.dime.pixdim(2:4) = voxel_size;
  hdr.hist.originator(1:3) = originator;

  %  set sform or qform to non-use, because they have been
  %  applied in xform_nii
  %
  hdr.hist.qform_code = 0;
  hdr.hist.sform_code = 0;
end

%  apply space_unit to pixdim if not 1 (mm)
%
space_unit = get_units(hdr);

if space_unit ~= 1
  hdr.dime.pixdim(2:4) = hdr.dime.pixdim(2:4) * space_unit;

  %  set space_unit of xyzt_units to millimeter, because
  %  voxel_size has been re-scaled
  %
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,1,0));
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,2,1));
  hdr.dime.xyzt_units = char(bitset(hdr.dime.xyzt_units,3,0));
end

return;         % change_hdr


%-----------------------------------------------------------------------
function orient = get_orient(R)

orient = [];

for i = 1:3
  switch find(R(i,:)) * sign(sum(R(i,:)))
    case 1
      orient = [orient 1]; %#ok   % Left to Right
    case 2
      orient = [orient 2]; %#ok   % Posterior to Anterior
    case 3
      orient = [orient 3]; %#ok   % Inferior to Superior
    case -1
      orient = [orient 4]; %#ok   % Right to Left
    case -2
      orient = [orient 5]; %#ok   % Anterior to Posterior
    case -3
      orient = [orient 6]; %#ok   % Superior to Inferior
  end
end

return;         % get_orient


%-----------------------------------------------------------------------
function [space_unit, time_unit] = get_units(hdr)

switch bitand(hdr.dime.xyzt_units, 7) % mask with 0x07
  case 1
    space_unit = 1e+3;    % meter, m
  case 3
    space_unit = 1e-3;    % micrometer, um
  otherwise
    space_unit = 1;     % millimeter, mm
end

switch bitand(hdr.dime.xyzt_units, 56)  % mask with 0x38
  case 16
    time_unit = 1e-3;     % millisecond, ms
  case 24
    time_unit = 1e-6;     % microsecond, us
  otherwise
    time_unit = 1;      % second, s
end

return;         % get_units

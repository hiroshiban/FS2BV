function ImportFreeSurfer2BrainVoyager(fs_subj_dir,do_flg)

% Imports anatomy, auto-segmentation, surface, and volume/surface ROI files generated by FreeSurfer into BrainVoyager
% function ImportFreeSurfer2BrainVoyager(fs_subj_dir,:do_flg)
%
% This function reads anatomy, auto-segmentation, surface, and volume/surface
% ROI files and converts them to BrainVoyager VMR, SRF, and VOI/POI data.
%
% For FreeSurfer details, please see
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki
% https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
% https://surfer.nmr.mgh.harvard.edu/fswiki/SubcorticalSegmentation
% https://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011
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
% This ImportFreeSurfer2BrainVoyager function especially uses
% 1. ConvertFreeSurferMGZ2VMR for importing anatomy and auto segmentation files as VMRs
% 2. ConvertFreeSurferSurf2SRF for importing surface files as SRFs
% 3. ConvertFreeSurferParcellation2BVvoi for importing volume ROI files as VOIs
% 4. ConvertFreeSurferAnnotation2BVpoi for importing surface ROI files as POIs
%
% [example]
% >> fs_subj_dir='./freesurfer/subjects/DC';
% >> do_flg=ones(1,5); % run all the importing procedures
% >> ImportFreeSurfer2BrainVoyager(fs_subj_dir,do_flg);
%
% [input]
% fs_subj_dir : FreeSurfer subject directory with a relative path format, in which
%               the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC
% do_flg      : (optional) a [1 x 5] matrix. flags to decide whether running the
%               importing step listed below.
%               step 01: importing FreeSurfer anatomy files
%               step 02: importing FreeSurfer volume ROI files
%               step 03: importing FreeSurfer surface files
%               step 04: importing FreeSurfer surface ROI files
%               step 05: adding subject name prefix (e.g. HB_) at the head of each of files.
%               if each flag is set to 1, the importing procedure correspond to that step
%               is run. If some of the flags are set to 0, the corresponding procure(s)
%               will be skipped.
%               do_flg=[1,1,1,1,1]; by default.
%
% [output]
% no output variable
% generated BrainVoyager-format files are stored in fs_subj_dir as
% fs_subj_dir/BrainVoyager/(subj_name)/{vmr|srf|voi|poi}/*.{vmr|srf|voi|poi}
%
% [dependency]
% 1. freesurfer_matlab_tools
% Matlab tools to read/write FreeSurfer files.
% The tools are copied under ~/freesurfer/matlab/ when you install FreeSurfer.
%
% 2. geom3d
% geom3d library is to handle and visualize 3D geometric primitives
% such as points, lines, planes, polyhedra... It provides low-level functions
% for manipulating 3D geometric primitives, making easier the development of more
% complex geometric algorithms.
%
% 3. UNIX command line tools such as gzip are required to use this function.
% Especially, if you want to use this function on Windows, please install UNIX command emulator,
% Cygwin or MSYS2, and set an envitonmental path to the tools in advance.
%
%
% Created    : "2017-09-11 18:40:01 ban"
% Last Update: "2021-06-16 10:01:52 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(do_flg), do_flg=[1,1,1,1,1]; end

% check the freesurfer directory and files
if ~exist(fullfile(pwd,fs_subj_dir),'dir')
  error('fs_subj_dir not found. check the input variable.');
end

% check whether the do_flg is set properly.
if numel(do_flg)~=5
  error('do_flg is a [1 x 5] matrix. please check the input variable.');
end

%% preparing a directory to store the converted files
[dummy,subj]=fileparts(fullfile(pwd,fs_subj_dir)); %#ok
save_dir=fullfile(pwd,fs_subj_dir,'BrainVoyager');
if ~exist(save_dir,'dir'), mkdir(save_dir); end
save_dir=fullfile(pwd,fs_subj_dir,'BrainVoyager',subj);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

dir_strings={'vmr','voi','srf','poi'};
for dd=1:1:length(dir_strings)
  if ~exist(fullfile(save_dir,dir_strings{dd}),'dir'), mkdir(fullfile(save_dir,dir_strings{dd})); end
end

%% set the paths to the tools for handling FreeSurfer files and manipulating 3D meshes on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));
addpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','geom3d'));
addpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','meshes3d'));

fprintf('\nImporting FreeSurfer files into BrainVoyager...\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 01: importing FreeSurfer anatomy files
%           usage: function vmr=ConvertFreeSurferMGZ2VMR(mgz_file,VMR_coregistration_file,ACPC_file,TAL_file,MNI_file,asis_conversion_flg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*************************************************************\n');
fprintf('step 01: importing FreeSurfer anatomy files...');
if do_flg(1)
  % importing T1.mgz
  mgz_file=relativepath(fullfile(pwd,fs_subj_dir,'mri','T1.mgz'));
  if strcmpi(mgz_file(end),filesep()), mgz_file=mgz_file(1:end-1); end % end-1 is required to omit the filesep() at the end of the string
  if ~exist(mgz_file,'file')
    fprintf('no T1.mgz file found. skipping...');
  else
    string=evalc('ConvertFreeSurferMGZ2VMR(mgz_file);');
    movefile([mgz_file,'.vmr'],fullfile(save_dir,dir_strings{1}));
  end

  % importing ribbon.mgz
  ribbon_file=relativepath(fullfile(pwd,fs_subj_dir,'mri','ribbon.mgz'));
  if strcmpi(ribbon_file(end),filesep()), ribbon_file=ribbon_file(1:end-1); end % end-1 is required to omit the filesep() at the end of the string
  if ~exist(ribbon_file,'file')
    fprintf('no ribbon file found. skipping...');
  else
    string2=evalc('ConvertFreeSurferRibbon2BL2VMR(ribbon_file);');
    movefile([ribbon_file,'_LH_BL2.vmr'],fullfile(save_dir,dir_strings{1}));
    movefile([ribbon_file,'_RH_BL2.vmr'],fullfile(save_dir,dir_strings{1}));
  end
  fprintf('completed.\n');
else
  fprintf('skipped.\n');
end
fprintf('*************************************************************\n');
fprintf('\n');
if do_flg(1), disp(string); disp(string2); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 02: importing FreeSurfer volume ROI files
%           usage: function voi=ConvertFreeSurferParcellation2BVvoi(fs_subj_dir,VMR_coregistration_file,ACPC_file,TAL_file,MNI_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*************************************************************\n');
fprintf('step 02: importing FreeSurfer volume ROI files...');
if do_flg(2)
  string=evalc('ConvertFreeSurferParcellation2BVvoi(fs_subj_dir);');
  movefile(fullfile(pwd,fs_subj_dir,'mri','*.voi'),fullfile(save_dir,dir_strings{2}));
  fprintf('completed.\n');
else
  fprintf('skipped.\n');
end
fprintf('*************************************************************\n');
fprintf('\n');
if do_flg(2), disp(string); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 03: importing FreeSurfer surface files
%           usage: function ConvertFreeSurferSurf2SRF(fs_subj_dir,VMR_coregistration_file,display_flg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*************************************************************\n');
fprintf('step 03: importing FreeSurfer surface files...');
if do_flg(3)
  string=evalc('ConvertFreeSurferSurf2SRF(fs_subj_dir,'''',0);');
  movefile(fullfile(pwd,fs_subj_dir,'surf','*.srf'),fullfile(save_dir,dir_strings{3}));
  % re-attaching reference SRFs
  string2=evalc(['hemis={''lh'',''rh''};',...
                 'for hh=1:1:length(hemis),',...
                 '  srffiles=GetFiles(fullfile(pwd,fs_subj_dir,''BrainVoyager'',subj,''srf''),'...
                 '                    ''*.srf'',{sprintf(''%s'',hemis{hh});{''orig'',''pial'',''smoothwm'',''white''}});',...
                 '  for ii=1:1:length(srffiles),',...
                 '    srf=BVQXfile(srffiles{ii});',...
                 '    srf.AutoLinkedMTC=fullfile(pwd,fs_subj_dir,''BrainVoyager'',subj,''srf'',sprintf(''%s_%s.smoothwm.srf'',subj,hemis{hh}));',...
                 '    srf.Save();',...
                 '    srf.ClearObject();',...
                 '    clear srf;',...
                 '  end;',...
                 'end;']);
  fprintf('completed.\n');
else
  fprintf('skipped.\n');
end
fprintf('*************************************************************\n');
fprintf('\n');
if do_flg(3), disp(string); disp(string2); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 04: importing FreeSurfer surface ROI files
%           usage: function poi=ConvertFreeSurferAnnotation2BVpoi(fs_subj_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*************************************************************\n');
fprintf('step 04: importing FreeSurfer surface ROI files...');
if do_flg(4)
  string=evalc('ConvertFreeSurferAnnotation2BVpoi(fs_subj_dir);');
  movefile(fullfile(pwd,fs_subj_dir,'label','*.annot.poi'),fullfile(save_dir,dir_strings{4}));
  fprintf('completed.\n');
else
  fprintf('skipped.\n');
end
fprintf('*************************************************************\n');
fprintf('\n');
if do_flg(4), disp(string); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 05: adding subject name at the head of each of files.
%           usage: function filenames=AddPrefix2Filename(target_dir,extension,:prefix_tgt,:prefix_head,:prefix_tail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n*************************************************************\n');
fprintf('step 05: adding subject name at the head of each of files...');
if do_flg(5)
  string=cell(5,1);
  tgt_str={'poi','srf','vmr','voi'};
  for ii=1:1:length(tgt_str)
    string{ii}=evalc('AddPrefix2Filename(fullfile(fs_subj_dir,''BrainVoyager'',subj,tgt_str{ii}),''*'',{'''';{sprintf(''%s_'',subj)}},sprintf(''%s_'',subj),'''');');
  end
  fprintf('completed.\n');
else
  fprintf('skipped.\n');
end
fprintf('*************************************************************\n');
fprintf('\n');
if do_flg(5)
  for ii=1:1:length(string), disp(string{ii}); end
end


% remove the paths to the tools for handling FreeSurfer files and manipulating 3D meshes on MATLAB
if ~isempty(strfind(path,'freesurfer_matlab_tools'))
  rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));
end
if ~isempty(strfind(path,['geom3d',filesep(),'geom3d']))
  rmpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','geom3d'));
end
if ~isempty(strfind(path,['geom3d',filesep(),'meshes3d']))
  rmpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','meshes3d'));
end

return

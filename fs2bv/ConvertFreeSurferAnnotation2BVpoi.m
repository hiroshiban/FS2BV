function poi=ConvertFreeSurferAnnotation2BVpoi(fs_subj_dir)

% Reads FreeSurfer surface annotation files (*.annot) and converts them to BrainVoyager POI files.
% function poi=ConvertFreeSurferAnnotation2BVpoi(fs_subj_dir)
%
% This function reads FreeSurfer annotation files, *.annot (those files are in ~/fs_subj_dir/label/),
% and converts them to BrainVoyager POI files.
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
% To find the VOIs in which specific XYZ coordinates in TAL/MNI space, please use the function below,
% GetAreaNameFromAtlasVOI               : Returns area candidates, in which the input XYZ coordinate(s)
%                                        is(are) belonging to, based on the pre-defined VOI atlases.
%
% [example]
% >> fs_subj_dir='./freesurfer/subjects/DC';
% >> poi=ConvertFreeSurferAnnotation2BVpoi(fs_subj_dir);
%
% [input]
% fs_subj_dir : FreeSurfer subject directory with a relative path format, in which
%               the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC
%
% [output]
% poi         : BrainVoyager POI object, a cell structure, poi{hemis}{num_annot_files}
% The generated POI files are saved in the same directory with the input auto annotation files
% with an extention '.poi'.
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
% Created    : "2017-09-11 13:40:01 ban"
% Last Update: "2018-09-04 18:47:33 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end

% check the freesurfer directory and files
if ~exist(fullfile(pwd,fs_subj_dir),'dir')
  error('fs_subj_dir not found. check the input variable.');
end

% set files to be processed
hemis={'lh','rh'};
for hh=1:1:length(hemis)
  fsfiles{hh}=GetFiles(fullfile(pwd,fs_subj_dir),sprintf('%s.*.annot',hemis{hh}));

  % omitting irrelevant files
  fscounter=0;
  fs_idx=logical(ones(length(fsfiles{hh}),1));
  for ii=1:1:length(fsfiles{hh})
    [dummy1,dummy2,fsext]=fileparts(fsfiles{hh}{ii});
    if strcmpi(fsext,'.poi')
      fs_idx(ii)=false;
    end
  end
  fsfiles{hh}=fsfiles{hh}(fs_idx);
end

% set the path to the tools for handling FreeSurfer files on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

% display message
fprintf('Target FreeSurfer subject directory : %s\n',fs_subj_dir);
fprintf('\n');

%% processing
for hh=1:1:length(hemis)
  poi{hh}=cell(length(fsfiles{hh}),1);
  for ii=1:1:length(fsfiles{hh})

    % display message
    [fspath,fsfname,fsext]=fileparts(fsfiles{hh}{ii});
    fprintf('loading a FreeSurfer annotation file onto MATLAB: %s%s...',fsfname,fsext);

    % loading surface annotation file
    [vertex,label,clut]=read_annotation(fsfiles{hh}{ii});

    if strcmpi(clut.struct_names{1},'Unknown')
      clut.numEntries=clut.numEntries-1;
      clut.struct_names=clut.struct_names(2:end);
      clut.table=clut.table(2:end,:);
    end

    fprintf('done.\n');

    % initializing a BrainVoyager POI file
    poi{hh}{ii}=BVQXfile('new:poi');
    poi{hh}{ii}.NrOfMeshVertices=numel(vertex);
    poi{hh}{ii}.NrOfPOIs=clut.numEntries;
    for pp=2:1:clut.numEntries, poi{hh}{ii}.POI(pp)=poi{hh}{ii}.POI(1); end % just initialization

    % converting them to a BrainvVoyager POI file
    for pp=1:1:clut.numEntries
      fprintf('converting to POI (%02d): %s...',pp,clut.struct_names{pp});
      poi{hh}{ii}.POI(pp).Name=clut.struct_names{pp};
      poi{hh}{ii}.POI(pp).Color=clut.table(pp,1:3);

      poi{hh}{ii}.POI(pp).Vertices=find(label==clut.table(pp,5));
      poi{hh}{ii}.POI(pp).NrOfVertices=numel(poi{hh}{ii}.POI(pp).Vertices);
      if poi{hh}{ii}.POI(pp).NrOfVertices>0
        poi{hh}{ii}.POI(pp).LabelVertex=poi{hh}{ii}.POI(pp).Vertices(round(poi{hh}{ii}.POI(pp).NrOfVertices/2));
      else
        poi{hh}{ii}.POI(pp).LabelVertex=0;
      end
      fprintf('done.\n');
    end

    fprintf('saving: %s...',[fsfname,fsext,'.poi']);
    poi{hh}{ii}.SaveAs(fullfile(fspath,[fsfname,fsext,'.poi']));
    fprintf('done.\n');
    fprintf('\n');

  end % for ii=1:1:length(fsfiles)
end % for hh=1:1:length(hemis)

% clean up
if ~nargout
  for hh=1:1:length(hemis)
    for ii=1:1:length(fsfiles{hh})
      if ~isempty(poi{hh}{ii})
        poi{hh}{ii}.ClearObject();
      end
    end
  end
  clear poi;
end

fprintf('completed.\n');

% remove the path to the tools for handling FreeSurfer files on MATLAB
rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));

return

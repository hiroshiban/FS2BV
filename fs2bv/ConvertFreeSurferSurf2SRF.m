function ConvertFreeSurferSurf2SRF(fs_subj_dir,VMR_coregistration_file,display_flg)

% Reads FreeSurfer surface files and converts them to BrainVoyager SRF files.
% function ConvertFreeSurferSurf2SRF(fs_subj_dir,:VMR_coregistration_file,:display_flg)
% (: is optional)
%
% This function reads FreeSurfer reconstructed cortical surface files, {lh|rh}.orig, {lh|rh}.white,
% {lh|rh}.smoothwm, {lh|rh}.pial, {lh|rh}.inflated, and {lh|rh}.sphere (those files are in ~/fs_subj_dir/surf/),
% and converts them to BrainVoyager SRF files.
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
% [example]
% >> fs_subj_dir='./freesurfer/subjects/DC';
% >> ConvertFreeSurferParcellation2BVvoi(fs_subj_dir);
%
% [input]
% fs_subj_dir : FreeSurfer subject directory with a relative path format, in which
%               the origin of the path is the location where this function is called.
%               e.g. ../subjects/DC
% VMR_coregistration_file : (optional) a transformation matrix file that is generated after coregistering
%               (fs_subj_dir)/mri/wm.seg.mgz(vmr) to the corresponding subject's BrainVoyager VMR
%               file in the native space (generally *_final.vmr in BVQX_hbtools).
%               for the details of coregistration, please see MaskVMRbyFreeSurferSegmentation.m
%               The file should be specified as a relative path format.
% display_flg : (optional) whether displaying the converted surface. [0|1].
%               if set to 1, the generated SRF surface is displayed on MATLAB figure window.
%               0 by default.
%
% [output]
% The generated SRF files are saved in the same directory with the input surface files
% with an extention '.srf'.
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
% Created    : "2017-09-11 09:00:48 ban"
% Last Update: "2021-06-16 10:02:39 ban"

%% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(VMR_coregistration_file), VMR_coregistration_file=''; end
if nargin<3 || isempty(display_flg), display_flg=0; end

%% check the freesurfer directory and files
if ~exist(fullfile(pwd,fs_subj_dir),'dir')
  error('fs_subj_dir not found. check the input variable.');
end

if ~isempty(VMR_coregistration_file) && ~exist(fullfile(pwd,VMR_coregistration_file),'file')
  error('VMR_coregistration_file not found. check the input variable.');
end

% set files to be processed
hemis={'lh','rh'};
fsfiles=cell(6,2); % 6 = orig, white, smoothwm, pial, inflated, and sphere, 2 = lh and rh;
cvfiles=cell(6,2); % 6 = orig, white, smoothwm, pial, inflated, and sphere, 2 = lh and rh
for hh=1:1:length(hemis)
  fsfiles{1,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.orig',hemis{hh})); % FreeSurfer Surface file
  fsfiles{2,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.white',hemis{hh}));
  fsfiles{3,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.smoothwm',hemis{hh}));
  fsfiles{4,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.pial',hemis{hh}));
  fsfiles{5,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.inflated',hemis{hh}));
  fsfiles{6,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.sphere',hemis{hh}));

  cvfiles{1,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv',hemis{hh})); % curvature file
  cvfiles{2,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv',hemis{hh}));
  cvfiles{3,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv',hemis{hh}));
  cvfiles{4,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv.pial',hemis{hh}));
  cvfiles{5,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv',hemis{hh}));
  cvfiles{6,hh}=fullfile(pwd,fs_subj_dir,'surf',sprintf('%s.curv',hemis{hh}));
end

for hh=1:1:length(hemis)
  for ii=1:1:length(fsfiles)
    if ~exist(fsfiles{ii,hh},'file')
      %[dummy,fsfname,fsext]=fileparts(fsfiles{ii,hh});
      %warning('%s%s not found. omitting it from the processing',fsfname,fsext);
      fsfiles{ii,hh}='';
    end

    if ~exist(cvfiles{ii,hh},'file')
      %[dummy,cvfname,cvext]=fileparts(cvfiles{ii,hh});
      %warning('%s%s not found. omitting it from the processing',cvfname,cvext);
      cvfiles{ii,hh}='';
    end
  end
end

% set the origin (mesh center) of the surface.
% Here, in FreeSurfer, while the surface files treat the 128,128,128 voxel
% (of the T1.mgz and other images) as the origin, the center of the surface
% is [0,0,0]. It means that the surface origin [0,0,0] seems to be registered
% to voxel [128,128,128] internally.
% Then, in this function, mesh center should be set to [0,0,0], not [128,128,128]
mcenter=[0,0,0];
mshift=[128,128,128]; % mshift is used to match the surface coordinates with the original VMR (e.g. T1.mgz.vmr) file.

% add paths to the tools for handling FreeSurfer files and manipulating 3D meshes on MATLAB
addpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));
addpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','geom3d'));
addpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','meshes3d'));

% display message
fprintf('Target FreeSurfer subject directory : %s\n',fs_subj_dir);
fprintf('\n');

%% processing
for hh=1:1:length(hemis)
  for ii=1:1:length(fsfiles)

    if isempty(fsfiles{ii,hh}), continue; end

    % display message
    [fspath,fsfname,fsext]=fileparts(fsfiles{ii,hh});
    fprintf('converting FreeSurfer surface to BrainVoayger SRF: %s%s --> %s%s.srf...',fsfname,fsext,fsfname,fsext);

    [vertex,faces]=read_surf(fsfiles{ii,hh});
    vertex=vertex(:,[1,3,2]); % triangular order should be flipped (rotated) between left/right as FS and BV takes the different coordinate system

    % affine transformation
    % Since FreeSurfer and BrainVoyager take different coordinate systems,
    % the vertex xyz coordinates should be rotated and flipped as below.

    % here, the rotation along y-axis corresponds to that along z-axis
    %       the rotation along z-axis corresponds to that along x-axis
    vertex=transformVector3d(vertex-repmat(mcenter,[size(vertex,1),1]),createRotationOy(pi/2))+repmat(mcenter,[size(vertex,1),1]);
    vertex=transformVector3d(vertex-repmat(mcenter,[size(vertex,1),1]),createRotationOz(pi))+repmat(mcenter,[size(vertex,1),1]);

    % shift vertex coordinates to match them with the VMR coordinates.
    vertex=vertex+repmat(mshift,[size(vertex,1),1]);

    % for the inflated surfaces, I will shift the vertex along x-axis,
    % that means z-shift in FreeSufer coordinates on BrainVoyager.
    % this is required as otherwise the left and right surfaces are
    % overlapped when they are loaded in the same window simultaneously.
    % note: this coordinate shift will never affect the mapping results
    % as the mapping on the inflated surfaces are done using the reference
    % surfaces' (the original white matter surfaces) coordinates.
    if strcmpi(fsext,'.inflated'), vertex=vertex+repmat([0,0,(-1)^(hh-1)*64],[size(vertex,1),1]); end

    % converting FreeSurfer surf to BrainVoyager SRF
    srf=BVQXfile('new:srf');
    srf.MeshCenter=mshift;
    srf.NrOfTriangles=size(faces,1);
    srf.NrOfVertices=size(vertex,1);
    srf.VertexCoordinate=vertex;
    srf.TriangleVertex=faces+1;
    srf.VertexNormal=vertexNormal(vertex,faces+1);
    warning off; %#ok
    try
      srf.Neighbors=srf.TrianglesToNeighbors();
    catch
      % Sometimes srf.TrianglesToNeighbors can not compute the neighbors properly due to
      % the inconsistency of the definition of the 'closed' surfaces between FreeSurfer and
      % BrainVoyager. Therefore, I prepared a modified version of mesh_trianglestoneighbors,
      % in which validation steps of the bad vertices and their neighbors are skipped to avoid
      % faling into an inifinite loop. Furthermore, the function is modified so that it stops
      % the loops even when the neighbors can not be found after 1000 repetitions for each of nodes.
      fprintf('\nWARNING: found inconsistency between the FreeSurfer and BrainVoyager surfaces. ignoring surface hole(s). please be careful.\n');
      srf.Neighbors=mesh_trianglestoneighborsNoCheck(srf.NrOfVertices,srf.TriangleVertex);
    end
    warning on;
    srf.VertexColor=repmat(zeros(1,4),srf.NrOfVertices,1); % initializing vertex colors with a reference coloring mode
    % setting curvature colors with a direct coloring mode
    if ~isempty(cvfiles{ii,hh})
      curv=read_curv(cvfiles{ii,hh});
      gyidx=find(curv<0);
      slidx=find(curv>=0);
      srf.VertexColor(gyidx,:)=repmat([NaN,175,175,175],[numel(gyidx),1]);
      srf.VertexColor(slidx,:)=repmat([NaN,100,100,100],[numel(slidx),1]);
      clear curv gyidx slidx;
    end
    srf.ConvexRGBA(4)=1;
    srf.ConcaveRGBA(4)=1;

    % attach the reference surface for inflated and sphere surfaces.
    if ii>=5, srf.AutoLinkedMTC=[fsfiles{3,hh},'.srf']; end

    % apply affine transformation if VMR_coregistration_file is given
    if ii<6 && ~isempty(VMR_coregistration_file) % if the target surface is not a morphed sphere.
      trf=BVQXfile(fullfile(pwd,VMR_coregistration_file));
      srf=srf.Transform(trf.TFMatrix);
      trf.ClearObject(); clear trf;
    end

    % saving
    srf.SaveAs(fullfile(fspath,sprintf('%s%s.srf',fsfname,fsext)));
    srf.ClearObject(); clear srf;

    fprintf('done.\n');
  end % for ii=1:1:length(fsfiles)
end % for hh=1:1:length(hemis)

% for visualization
if display_flg
  for hh=1:1:length(hemis)
    for ii=1:1:length(fsfiles)
      if isempty(fsfiles{ii,hh}), continue; end
      srffile=relativepath([fsfiles{ii,hh},'.srf']);
      if strcmpi(srffile(end),filesep()), srffile=srffile(1:end-1); end
      SurfViewerSRF(srffile,[],0);
    end
  end
end

% remove paths to the tools for handling FreeSurfer files and manipulating 3D meshes on MATLAB
rmpath(fullfile(fileparts(mfilename('fullpath')),'freesurfer_matlab_tools'));
rmpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','geom3d'));
rmpath(fullfile(fileparts(mfilename('fullpath')),'geom3d','meshes3d'));

return


%%% subfunction

% The sub-function below is modified by H.Ban for compatibilities with BrainVoyager
% file format from a function distributed by Dr Don Hagler with his courtesy,

function nbrs = fs_find_neighbors(faces,nverts)

% nbrs = fs_find_neighbors(surf)
%
% Inputs:
%   faces:  vertex numbers for each face (3 corners)
%   nverts: number of vertices contained in the input cortex
%
% Output:
% surf is a structure containg:
%   nbrs:   vertex numbers of neighbors for each vertex {num_vertex,2}
%           nbrs{n,1}=num_vertices
%           nbrs{n,2}=vertices IDs
%
% original code for finding neighbors taken from Moo Chung's mni_getmesh
% modified from fs_find_neighbors, written by Dr Don Hagler 05/09/06
%
%
% Created    : "ban 2018-05-18 16:52:50"
% Last Update: "ban 2018-05-18 16:52:50"

nfaces=size(faces,1);

% compute the maximum degree of node -- number of edges = number of neighbors
num_nbrs=zeros(nverts,1);
for ii=1:1:nfaces, num_nbrs(faces(ii,:))=num_nbrs(faces(ii,:))+1; end
max_num_nbrs=max(num_nbrs);

% find nearest neighbors
ns=zeros(nverts,max_num_nbrs);
nbrs=cell(nverts,2);
for vv=1:1:nverts, nbrs{vv,1}=0; end
for ii=1:1:nfaces
  for jj=1:1:3
    vcur=faces(ii,jj);
    if jj==1
      vidx=[2,3];
    elseif jj==2
      vidx=[3,1];
    else
      vidx=[1,2];
    end
    for vv=vidx%1:3
      vnbr=faces(ii,vv);
      if find(ns(vcur,:)==vnbr), continue; end
      ns(vcur,min(find(ns(vcur,:)==0)))=vnbr;
      nbrs{vcur,1}=nbrs{vcur,1}+1;
      nbrs{vcur,2}(nbrs{vcur,1})=vnbr;
      %if numel(nbrs{1,2})<5, nbrs{1,2},jj, end
    end
  end
end
%keyboard

return

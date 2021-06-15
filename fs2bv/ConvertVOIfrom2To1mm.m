function voi=ConvertVOIfrom2To1mm(VOI_dir,prefix_voi,save_flg)

% Converts VOI(s) defined in 2.0 mm space into standard 1.0 mm VOI(s).
% function voi=ConvertVOIfrom1To0p5mm(VOI_dir,:prefix_voi,:save_flg)
% (: is optional)
%
% This function reads VOI(s) defined in 2.0 mm cubic space (256 dimensions) and
% converts it(them) into BrainVoyager default 1.0 mm VOI(s) (256 dimensions) by
% upsampling the neighboring voxels. This upsampling will be useful in importing
% SPM, FSL, or FreeSurfer ROIs into BrainVoyager since some ROIs used for these
% software package are defined in 2.0 mm, not 1.0 mm, cubic space.
%
% [input]
% VOI_dir    : Target directory that contains 1.0 mm VOI files
%              e.g. '\CD\zk08_382\ROI_vtc\voi_files'
%              Target directory should be specified as such
%              the current directory where this function is
%              called is the origin
% prefix_voi : (optional) string to determine the target VOI
%              from multiple files, e.g. prefix_vmp='all_combined';
%              empty by default.
% save_flg   : (optional) whther saving the generated VOI. 1 by default.
%              if 1, the generated VOI will be saved with "_highres" prefix.
%
% [output]
% voi        : BrainVoyager VOI object(s) defined in 1.0 mm (256 dimensions) space
%              a cell structure
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
%
% Created    : "2018-08-27 11:18:21 ban"
% Last Update: "2018-08-27 11:53:18 ban"

% check the input variables
if nargin<1, help(mfilename()); return; end
if nargin<2 || isempty(prefix_voi), prefix_voi=''; end
if nargin<3 || isempty(save_flg), save_flg=1; end

% check the directory
if ~exist(fullfile(pwd,VOI_dir),'dir'), error('can not find VOI_dir. check the input variable'); end
fprintf('TARGET: %s\n',fullfile(pwd,VOI_dir));

% get the VOI files
voifiles=GetFiles(fullfile(pwd,VOI_dir),'*.voi',prefix_voi);
if isempty(voifiles), error('can not find any VOI file in VOI_dir. check the input variables.'); end

% processing
if nargout, voi=cell(length(voifiles),1); end
for ii=1:1:length(voifiles)
  [voipath,voifname,voiext]=fileparts(voifiles{ii});
  fprintf('processing: %s%s...',voifname,voiext);

  svoi=BVQXfile(voifiles{ii});

  if svoi.OriginalVMRResolutionX~=1 || svoi.OriginalVMRResolutionY~=1 || svoi.OriginalVMRResolutionZ~=1
    fprintf('warning: the VOI file is not defined in 1.0 mm space...skipping...\n');
    continue
  end

  hvoi=BVQXfile('new:voi');
  hvoi.FileVersion=4;
  hvoi.ReferenceSpace=svoi.ReferenceSpace;
  hvoi.OriginalVMROffsetX=0;
  hvoi.OriginalVMROffsetY=0;
  hvoi.OriginalVMROffsetZ=0;
  hvoi.OriginalVMRFramingCubeDim=256;
  hvoi.SubjectVOINamingConvention=svoi.SubjectVOINamingConvention;
  hvoi.Convention=svoi.Convention;
  hvoi.NrOfVOIs=svoi.NrOfVOIs;
  hvoi.NrOfVTCs=svoi.NrOfVTCs;
  hvoi.VTCList=svoi.VTCList;
  hvoi.RunTimeVars=svoi.RunTimeVars;

  % converting 2mm voxels to 1mm ones by upsampling
  hvoi.VOI=svoi.VOI;
  for vv=1:1:svoi.NrOfVOIs
    voxels=zeros(8*size(svoi.VOI(vv).NrOfVoxels,1),3);
    basecoords=hvoi.VOI(vv).Voxels;
    for mm=1:1:size(basecoords,1)
      voxels(8*(mm-1)+1,:)=basecoords(mm,:)+[ 0, 0, 0];
      voxels(8*(mm-1)+2,:)=basecoords(mm,:)+[-1, 0, 0];
      voxels(8*(mm-1)+3,:)=basecoords(mm,:)+[-1,-1, 0];
      voxels(8*(mm-1)+4,:)=basecoords(mm,:)+[ 0,-1, 0];
      voxels(8*(mm-1)+5,:)=basecoords(mm,:)+[ 0, 0,-1];
      voxels(8*(mm-1)+6,:)=basecoords(mm,:)+[-1, 0,-1];
      voxels(8*(mm-1)+7,:)=basecoords(mm,:)+[ 0,-1,-1];
      voxels(8*(mm-1)+8,:)=basecoords(mm,:)+[-1,-1,-1];
    end
    % omit duplicated voxels
    hvoi.VOI(vv).Voxels=unique(voxels,'rows');
    hvoi.VOI(vv).NrOfVoxels=size(hvoi.VOI(vv).Voxels,1);
  end

  if save_flg
    hvoi.SaveAs(fullfile(voipath,[voifname,'_highres',voiext]));
  end

  if nargout, voi{ii}=hvoi; end

  fprintf('completed.\n');
end

return

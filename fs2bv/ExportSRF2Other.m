function ExportSRF2Other(SRF_dir,format_export,prefix_srf)

% Exports BVQX SRF file(s) to the other image fomat (vtk,byu,jvx).
% function ExportSRF2Other(SRF_dir,:format_export,:prefix_srf)
% (: is optional)
%
% [about]
% This function exports *.srf files to the other 3D graphics file format.
% supported formats are: VTK, BYU, JVX
%
% [exapmle]
% >> ExportSRF2Other('/HB/3D/zk09_034','VTK','LH')
%
% [input]
% SRF_dir       : target directory that contains SRF file (*.srf)
%                 e.g. '/HB/3D/zk09_034/
%                 the location should be specified as such
%                 the current directory where this function is
%                 called is the origin.
% format_export : (optional) file format you want to export.
%                 'VTK'(default), 'BYU', or 'JVX'
% prefix_srf    : (optoinal) string to determine the target srf
%                 from multiple files, e.g. 'HB'.
%                 empty by default.
%
% [output]
% no output variables
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
% [dependency]
% 1. BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
% 2. wildcardsearch.m
% enable reg-exp search of files
%
%
% Created    : "2010-06-10 21:22:49 ban"
% Last Update: "2018-09-03 10:59:44 ban"

% check input variables
if nargin<1, help(mfilename()); return ; end
if nargin<2 || isempty(format_export), format_export='VTK'; end
if nargin<3 || isempty(prefix_srf), prefix_srf=''; end

% check exporting file format
format_export=upper(format_export);
if ~strcmp(format_export,'VTK') && ~strcmp(format_export,'BYU') && ~strcmp(format_export,'JVX')
  error('Supported export format is ''VTK'', ''BYU'', or ''JVX''. Check format_export variable.');
end

% get srf files from the directories in the SRF_dir
srffiles=GetFiles(fullfile(pwd,SRF_dir),'*.srf',prefix_srf);

% display message
fprintf('Target SRF   : %s\n',fullfile(pwd,SRF_dir));

% attach a reference recosm file to each target srf file
for ii=1:1:length(srffiles)

  % display message
  [path,srfname,ext]=fileparts(srffiles{ii});
  fprintf('processing : %s%s --> %s.%s\n',srfname,ext,srfname,lower(format_export));

  % read srf file & create its object
  srf=BVQXfile(srffiles{ii});

  % export it as the other 3D graphics file fomat
  if ~strcmp(format_export,'VTK')
    eval(sprintf('srf.SaveAs%s(''%s.%s'');',format_export,[path,filesep(),srfname],lower(format_export)));
  else
    % export data as VTK format
    % !NOTE!
    % the default srf.SaveAsVTK in BVQXtools is insufficient to handle color-lookuptables and normals.
    % I thus wrote my own routines to convert SRF files to VTK files.

    vtkfile=[path,filesep(),srfname,'.vtk'];
    fid = fopen(vtkfile, 'wb', 'ieee-be');
    if fid < 1
      error( ...
          'BVQXfile:FileNotWritable', ...
          'Given filename is not writable.' ...
          );
    end

    % write header
    fwrite(fid, ['# vtk DataFile Version 3.0' char(10)], 'char');
    fwrite(fid, ['vtk output' char(10)], 'char');
    fwrite(fid, ['BINARY' char(10)], 'char');
    fwrite(fid, ['DATASET POLYDATA' char(10)], 'char');

    % write vertices
    numvert = srf.NrOfVertices;
    fwrite(fid, [sprintf('POINTS %d float', numvert) char(10)], 'char');
    fwrite(fid, (srf.VertexCoordinate' - srf.MeshCenter' * ones(1, numvert)), 'single');
    fwrite(fid, char(10), 'char');

    % write triangles
    numtri = srf.NrOfTriangles;
    fwrite(fid, [sprintf('POLYGONS %d %d', numtri, 4 * numtri) char(10)], 'char');
    fwrite(fid, [3 * ones(1, numtri); srf.TriangleVertex' - 1], 'int32');
    fwrite(fid, char(10), 'char');

    % write footer
    %fwrite(fid, [sprintf('CELL_DATA %d', numtri) char(10)], 'char');
    fwrite(fid, [sprintf('POINT_DATA %d', numvert) char(10)], 'char');
    fwrite(fid, char(10), 'char');

    % write vertex color
    colors=srf.VertexColor;
    colors=colors(:,[2:4,1]);
    colors(:,4)=255; % alpha

    fwrite(fid, [sprintf('COLOR_SCALARS scalars %d',4) char(10)], 'char');
    fwrite(fid, colors', 'uint8');
    fwrite(fid, char(10), 'char');

    % write normals
    fwrite(fid, ['NORMALS normals float' char(10)], 'char');
    fwrite(fid, srf.VertexNormal', 'single');
    fwrite(fid, char(10), 'char');

    % close file
    fclose(fid);
  end % if ~strcmp(format_export,'VTK')

  % clear BVQX object
  srf.ClearObject(); clear srf;

end % for ii=1:1:length(srffiles)

return

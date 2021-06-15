function ErasePOIBordersSRF(srf_as_input,srf_as_output,POI_file,POI_to_use,POI_not_to_use)

% Erases POI borders drawn on the SRF fle (regions to be erased are specified by the input POI file).
% function EraseCutBordersSRF(srf_as_input,srf_as_output,POI_file,:POI_to_use,:POI_not_to_use)
% (: is optional)
%
% This fuction removes POI borders drawn on the SRF file. The removal region(s)
% can be specified by the borders of POI_to_use and POI_file.
% This function can remove the borders more precisely than EraseCutBordersSRF
% and useful when you can delete exactly the borders of specific POI(s).
%
% !!! NOTICE !!!
% 1. The srf_as_input should be converted to direct color surface (gray-scale color one)
%    in advance of running this function.
% 2. The deleted borders are painted by white color. If you want to remove these white
%    lines, run Meshes -> Background and Curvature Colors... -> Curvature/Smooth
%
% [example]
% >> SRF_input='/HB/3d/zk09_034_LH_inflated_cut.srf';
% >> SRF_output='/HB/3d/zk09_034_LH_inflated_cut_2.srf';
% >> POI_file-='/HB/zk09_034/poi_files/all_ROIs_LH.poi';
% >> POI_to_use={'V1_lh','V2_lh'};
% >> POI_not_to_use={'V3d_lh','V3v_lh'};
% >> ErasePOIBordersSRF(SRF_input,SRF_output,POI_file,POI_to_use,POI_not_to_use);
%
% [input]
% srf_as_input   : target SRF file name with relative path
%                  e.g. '/HB/3D/zk09_034/zk09_034_LH_inflated.srf'
%                  the location should be specified as such
%                  the current directory where this function is
%                  called is the origin.
% srf_as_output  : output (border removed) SRF file name with relative path
%                  e.g. '/HB/3D/zk09_034/zk09_034_LH_inflated.srf'
%                  the location should be specified as such
%                  the current directory where this function is
%                  called is the origin.
% POI_file       : target SRF file name with relative path
%                  e.g. '/HB/zk09_034/poi_files/all_POI_lh.poi'
% POI_to_use     : (optional) POIs whose borders you want to cut from the
%                  target SRF. e.g. POI_to_use={'V1','V2'};
%                  if empty, all the POIs in POI_file will be used.
%                  empty by default
% POI_not_to_use : (optional) POIs whose borders should be prevent from
%                  being erased. the final borders to be erased are the borders
%                  that is included in POI_to_use but not included POI_not_to_use.
%                  For example, if you set
%                  POI_to_use={'V1','V2'};
%                  POI_not_to_use={'V3d','V3v'};
%                  then, the borders between V1 and V2 will be erased, whereas
%                  the borders between V2d and V3d, V2v and V3v will be preserved.
%                  if empty, no POI in POI_file will be used. empty by default
%
% [output]
% no output variable, the new border-removed SRF file will be saved as srf_as_output
%
% [dependency]
% 1. BVQXtools
% %MATLABtool%\BVQXtools_v08d is recommended.
%
%
% Created    : "2013-04-09 03:16:21 banh"
% Last Update: "2018-09-03 10:59:03 ban"

% check input variables
if nargin<3, help(mfilename()); return; end
if nargin<4 || isempty(POI_to_use), POI_to_use={}; end
if nargin<5 || isempty(POI_not_to_use), POI_not_to_use={}; end

% load the target SRF file
srffile=fullfile(pwd,srf_as_input);
if ~exist(srffile,'file'), error('srf_as_input is not found. check the input variable.'); end

% display message
[inputpath,inputfname,inputext]=fileparts(srf_as_input);
fprintf('Target SRF : %s%s\n',inputfname,inputext);
inputsrf=BVQXfile(srf_as_input);

% load the target POI file
poifile=fullfile(pwd,POI_file);
if ~exist(poifile,'file'), error('POI_file is not found. check the input variable.'); end

% display message
[poipath,poifname,poiext]=fileparts(poifile);
fprintf('Target POI : %s%s\n',poifname,poiext);
poi=BVQXfile(poifile);

% check whether the POIs are generated from the target SRF (simple version)
if inputsrf.NrOfVertices~=poi.NrOfMeshVertices
  error('vertices of the input SRF and POI mismatched. check the input variable.');
end

% check POI to be used
poicounter=0;
poilist=[];
if ~isempty(POI_to_use)
  for pp1=1:1:length(POI_to_use)
    for pp2=1:1:poi.NrOfPOIs
      if strcmp(poi.POI(pp2).Name,POI_to_use{pp1})
        poicounter=poicounter+1;
        poilist(poicounter)=pp2; %#ok
      end
    end
  end
  if length(POI_to_use)~=numel(poilist)
    disp('!!WARNING!! some POI(s) you specified in POI_to_use not found. check the input variable.');
  end
else
  poilist=1:poi.NrOfPOIs;
end

% check POI not to be used
poicounter=0;
notpoilist=[];
if ~isempty(POI_not_to_use)
  for pp1=1:1:length(POI_not_to_use)
    for pp2=1:1:poi.NrOfPOIs
      if strcmp(poi.POI(pp2).Name,POI_not_to_use{pp1})
        poicounter=poicounter+1;
        notpoilist(poicounter)=pp2; %#ok
      end
    end
  end
  if length(POI_not_to_use)~=numel(notpoilist)
    disp('!!WARNING!! some POI(s) you specified in POI_not_to_use not found. check the input variable.');
  end
else
  notpoilist=[];
end

% deleting borders on the surface

% get vertices to be erased
fprintf('\ndeleting borders...\n');
erase_vertices=[];
for pppp=poilist

  % display message
  fprintf('erasing borders   : %8s --> ',poi.POI(pppp).Name);

  % get SRF index within the target POI
  poivertices=(poi.POI(pppp).Vertices)';

  % get neighbors of the poivertices (nearest neighbor search).
  % this search is required as POI are defined as the inner regions of
  % cut borders. therefore, we need to search border vertices at neighboring
  % regions defined by vertices of POI one-by-one.
  nnpoivertices=[];
  for nn=poivertices
    nnpoivertices=[nnpoivertices,inputsrf.Neighbors{nn,2}]; %#ok
    nnpoivertices=unique(nnpoivertices);
  end

  % get border vertices
  bidx=find(inputsrf.VertexColor(nnpoivertices,2)==0 & ...
            inputsrf.VertexColor(nnpoivertices,3)==0 & ...
            inputsrf.VertexColor(nnpoivertices,4)==0);

  fprintf('% 4d vertices were found on the border(s).\n',numel(nnpoivertices(bidx)));

  erase_vertices=[erase_vertices,nnpoivertices(bidx)]; %#ok
end % for pppp=poilist
erase_vertices=unique(erase_vertices);

% get vertices to be preserved
preserve_vertices=[];
for pppp=notpoilist

  % display message
  fprintf('preserving borders: %8s --> ',poi.POI(pppp).Name);

  % get SRF index within the target POI
  poivertices=(poi.POI(pppp).Vertices)';

  % get neighbors of the poivertices (nearest neighbor search).
  % this search is required as POI are defined as the inner regions of
  % cut borders. therefore, we need to search border vertices at neighboring
  % regions defined by vertices of POI one-by-one.
  nnpoivertices=[];
  for nn=poivertices
    nnpoivertices=[nnpoivertices,inputsrf.Neighbors{nn,2}]; %#ok
    nnpoivertices=unique(nnpoivertices);
  end

  % get border vertices
  bidx=find(inputsrf.VertexColor(nnpoivertices,2)==0 & ...
            inputsrf.VertexColor(nnpoivertices,3)==0 & ...
            inputsrf.VertexColor(nnpoivertices,4)==0);

  fprintf('% 4d vertices were found on the border(s).\n',numel(nnpoivertices(bidx)));

  preserve_vertices=[preserve_vertices,nnpoivertices(bidx)]; %#ok
end % for pppp=notpoilist
preserve_vertices=unique(preserve_vertices);

% get the final vertices on the POI border(s) to be erased
erase_vertices=setdiff(erase_vertices,preserve_vertices);

% set no border flag on the surface
inputsrf.VertexColor(erase_vertices,1)=NaN;
inputsrf.VertexColor(erase_vertices,2:4)=255;
disp('completed deleting borders.');

% save the results
fprintf('saving the result as: %s...',srf_as_output);
inputsrf.SaveAs(fullfile(pwd,srf_as_output));
fprintf('done.\n');

% clean up
inputsrf.ClearObject(); clear inputsrf;
poi.ClearObject(); clear poi;

return

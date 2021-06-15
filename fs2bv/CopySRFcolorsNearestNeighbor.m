function CopySRFcolorsNearestNeighbor(srf_as_input,srf_to_copy,srf_as_output,RGB_to_copy)

% Copies RGB colors of one SRF file to the other one with a different structure by a Nearest-Neighbor sampling.
% function CopySRFcolorsNearestNeighbor(srf_as_input,srf_to_copy,srf_as_output,RGB_to_copy)
%
% This function copies the specified RGB color in one SRF(srf_as_input) surface to the other SRF(srf_to_copy)
% whose structure is different from the original one. Nearest-Neighbor samplings are used in copying.
% The generated SRF will be saved as srf_as_output.
%
% This function will be useful when you need to copy the ROI borders from
% one SRF to the other SRFs with some activations. To do that, please
% set RGB_to_copy as [0,0,0]
%
% [!!! NOTICE !!!] -- obsolete, Oct 31 2011, any SRF file can be acceptable now
% input SRF should be converted to direct colored surfaces in advance
% if srf_to_copy to which you want to copy the borders is the flat surface
%
% [input]
% srf_as_input  : input SRF file name 1 with relative path,
%                 should be specified so that the location where this
%                 function is called is the origin, e.g. './zk11_110/HB/zk11_110.polar.inflated_LH.srf'
% srf_to_copy   : input SRF file name 2 with relative path,
%                 should be specified so that the location where this
%                 function is called is the origin, e.g. './zk11_110/HB/zk11_110.flat_LH.srf'
% srf_as_output : output SRF file name with relative path,
%                 should be specified so that the location where this
%                 function is called is the origin, e.g. './zk11_110/HB/zk11_110.polar.flat_LH.srf'
% RGB_to_copy   : SRF's vertex color to copy, [R,G,B]
%
% [output]
% no output variable, the new SRF file with the specified color
% will be saved as srf_as_output
%
%
% Created    : "2014-06-21 17:13:40 ban"
% Last Update: "2018-09-04 18:47:54 ban"

% check input variables
if nargin<4, help(mfilename()); return; end
if numel(RGB_to_copy)~=3
  error('RGB_to_copy should be a matrix = [R,G,B]. check input variable');
end

% check whether two input SRFs exist
fprintf('checking input SRFs...');
fname_input=fullfile(pwd,srf_as_input);
fname_copy=fullfile(pwd,srf_to_copy);
fname_output=fullfile(pwd,srf_as_output);

if ~exist(fname_input,'file')
  error('srf_as_input is not found. check input variables.');
end

if ~exist(fname_copy,'file')
  error('srf_to_copy is not found. check input variables.');
end
fprintf('done.\n');

% generate BVQXobject
fprintf('generating BVQX SRF objects...');
input_srf=BVQXfile(fname_input);
copy_srf=BVQXfile(fname_copy);
fprintf('done.\n');

% get XYZ coordinates of the input/to_be_copied SRFs
if ~isempty(input_srf.AutoLinkedMTC)
   tmp=BVQXfile(input_srf.AutoLinkedMTC);
   input_xyz=tmp.VertexCoordinate;
   tmp.ClearObject(); clear tmp;
else
   input_xyz=input_srf.VertexCoordinate;
end

if ~isempty(copy_srf.AutoLinkedMTC)
   tmp=BVQXfile(copy_srf.AutoLinkedMTC);
   output_xyz=tmp.VertexCoordinate;
   tmp.ClearObject(); clear tmp;
else
   output_xyz=copy_srf.VertexCoordinate;
end

% copy the whole SRF structure to output_srf
fprintf('copying colors by nearest-neighbor samplings...');
output_srf=copy_srf.CopyObject();
copy_srf.ClearObject(); clear copy_srf;

% get the vertex IDs with the specified RGB value
%
% !!! NOTE !!!
% srf.VertexColor(:,1) holds vertex properties as below.
% 0, 1        : concave, convex
% 10001-10009 : cut/flattening regions
% 4.29..e+9   : empty(hidden) nodes for flattened surface
% NaN         : direct color map

%vertex_idx=find( input_srf.VertexColor(:,2)==RGB_to_copy(1) & ...
%                 input_srf.VertexColor(:,3)==RGB_to_copy(2) & ...
%                 input_srf.VertexColor(:,4)==RGB_to_copy(3) & ...
%                 (output_srf.VertexColor(:,1)<=20000 | isnan(output_srf.VertexColor(:,1))) );
vertex_idx=find( input_srf.VertexColor(:,2)==RGB_to_copy(1) & ...
                 input_srf.VertexColor(:,3)==RGB_to_copy(2) & ...
                 input_srf.VertexColor(:,4)==RGB_to_copy(3) );
vertex_idx=unique( ind2sub(size(input_srf.VertexColor),vertex_idx) ); % get rows

% copy the color of the vertices by Nearest-Neighbor samplings.
for vv=1:1:numel(vertex_idx)
  [dummy,vertex_nn_idx]=min( sum( ( output_xyz-repmat(input_xyz(vertex_idx(vv),:),[size(output_xyz,1),1]) ).^2, 2 ) );
  output_srf.VertexColor(vertex_nn_idx,:)=input_srf.VertexColor(vertex_idx(vv),:);
end

fprintf('done.\n');

% save it with srf_as_output name
fprintf('saving...');
output_srf.SaveAs(fname_output);
fprintf('done.\n');

% clean up
input_srf.ClearObject(); clear input_srf;
output_srf.ClearObject(); clear output_srf;

return

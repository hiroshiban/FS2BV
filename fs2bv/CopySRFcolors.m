function CopySRFcolors(srf_as_input,srf_to_copy,srf_as_output,RGB_to_copy)

% Copies RGB colors of one SRF file to the other (e.g. copy VOI boundary from inflated to flat).
% function CopySRFcolors(srf_as_input,srf_to_copy,srf_as_output,RGB_to_copy)
%
% This function copies the specified RGB color in one SRF(srf_as_input)
% surface to the other SRF(srf_to_copy) and save it as srf_as_output.
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
% Created    : "2011-07-13 09:34:17 banh"
% Last Update: "2018-09-04 18:47:51 ban"

% check input variables
if nargin<4, help(mfilename()); return; end
if numel(RGB_to_copy)~=3
  error('RGB_to_copy should be a matrix = [R,G,B]. check the input variable');
end

% check whether two input SRFs exist
fprintf('checking the input SRFs...');
fname_input=fullfile(pwd,srf_as_input);
fname_copy=fullfile(pwd,srf_to_copy);
fname_output=fullfile(pwd,srf_as_output);

if ~exist(fname_input,'file')
  error('srf_as_input is not found. check the input variables.');
end

if ~exist(fname_copy,'file')
  error('srf_to_copy is not found. check the input variables.');
end
fprintf('done.\n');

% generate BVQXobject
fprintf('generating BVQX SRF objects...');
input_srf=BVQXfile(fname_input);
copy_srf=BVQXfile(fname_copy);
fprintf('done.\n');

if input_srf.NrOfVertices~=copy_srf.NrOfVertices
  error('the number of vertices of the srf_as_input and srf_to_copy mismatched. check the input variables.');
end

% copy the whole SRF structure to output_srf
fprintf('copying colors...');
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

vertex_idx=logical( input_srf.VertexColor(:,2)==RGB_to_copy(1) & ...
                    input_srf.VertexColor(:,3)==RGB_to_copy(2) & ...
                    input_srf.VertexColor(:,4)==RGB_to_copy(3) & ...
                    (output_srf.VertexColor(:,1)<=20000 | isnan(output_srf.VertexColor(:,1))) );

% copy the color and save it with srf_as_output name
output_srf.VertexColor(vertex_idx,:)=input_srf.VertexColor(vertex_idx,:);

fprintf('done.\n');

fprintf('saving...');
output_srf.SaveAs(fname_output);
fprintf('done.\n');

% clean up
input_srf.ClearObject(); clear input_srf;
output_srf.ClearObject(); clear output_srf;

return

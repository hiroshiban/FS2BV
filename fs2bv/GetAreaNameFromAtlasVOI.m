function roi_name=GetAreaNameFromAtlasVOI(XYZ_coords,space,atlas_to_be_used)

% Returns ROI candidates, in which the input XYZ coordinate(s) is(are) belonging to, based on the pre-defined VOI atlases.
% function roi_name=GetAreaNameFromAtlasVOI(XYZ_coords,:space,:atlas_to_be_used)
% (: is optional)
%
% This function reads pre-defined VOI coordinates stored in ~/BVQX_hbtools/VOIs/,
% explores which VOIs the input XYZ coordinates are located, and returns the
% candidate ROI names. The output of this function will be useful to decide the
% VOI files to be used for the latter detailed ROI-based analyses.
%
% *TODO*
% to write a faster processing version by storing the VOI coods in a combined *.mat file.
%
% [example]
% >> roi_name=GetAreaNameFromAtlasVOI([66 -17 35],'MNI',[0,1,1,0,0,0,0,0,0]);
%
% [input]
% XYZ_coords  : XYZ coordinate(s) to be explored,  a [N (number of coordinates) x 3 (X,Y,Z)] matrix.
% space       : (optional) stereotaxic coordinate space, 'TAL' (Talairach) or 'MNI'.
%               'TAL' by default.
% atlas_to_be_used : (optional) atlas IDs to be used, [1 x 9] matrix consisted 0 and 1.
%               atlas_to_be_used=[1,1,1,1,1,1,1,1,1]; by default.
%               the atlases listed below with the corresponding IDs are set to 1 are
%               used for area search. the atlases whose IDs are 0 are omitted.
%               1. automated_anatomical_labeling_VOIs
%               2. Brainnetome_Atlas_VOIs
%               3. brodmann_VOIs
%               4. Buckner_JNeurophysiol11_MNI152_VOIs
%               5. Choi_JNeurophysiol12_MNI152_VOIs
%               6. FSL_atlases_VOIs
%               7. ProbAtlas_v4_VOIs
%               8. shenparcel_VOIs
%               9. Talairach.org_VOIs (when space='TAL' only)
%               when space='TAL', atlas_to_be_used(9) is always set to 0 since
%               'Talairach.org_VOIs' is defined only in TAL space.
%
% [output]
% roi_name    : ROI candidates found,
%               if multiple candidates are found, all the candidats are returned
%               as a cell structure consisted of the elements below.
%               roi_name{n}.name      : are name, e.g. V1
%                          .XYZ       : XYZ coordinate contained in the VOI.
%                          .atlas_type: atlas type, defined below as the "atlas_locations" variable.
%                          .voi_file  : the original VOI file that contains a VOI
%                                       in which XYZ_coords are included.
%
% [note]
% The area(s) of the input XYZ coordinates is(are) explored from the pre-defined VOI atlas stored in
% ~/BVQX_hbtools/VOIs/. If you want to add your own atlas (VOI files), please add the files into it
% and modify the "atlas_to_be_used" and "atlas_locations" variable-related codes below.
%
% [dependency]
% BVQXtools_v08d: BrainVoyager MATLAB tools.
%
%
% Created    : "2017-10-04 10:58:17 ban"
% Last Update: "2017-10-04 13:26:50 ban"

%% check the input variables
if nargin<1 || isempty(XYZ_coords), help(mfilename()); end
if nargin<2 || isempty(space), space='TAL'; end
if nargin<3 || isempty(atlas_to_be_used), atlas_to_be_used=ones(1,9); end

if size(XYZ_coords,2)~=3
  error('XYZ_coords should be [X_coord, Y_coord, Z_coord]. check the input variable.');
end

if ~strcmpi(space,'TAL') && ~strcmpi(space,'MNI')
  error('currently only ''TAL'' or ''MNI'' can be set to ''space''. check the input variable.');
end

if numel(atlas_to_be_used)~=9
  error('atlas_to_be_used should be a [1 x 9] matrix consisted of 0 and 1. check the input variable.');
end

if isempty(find(atlas_to_be_used~=0,1))
  error('to run the function properly, at least one of the elements of atlas_to_be_used should be 1. check the input variable.');
end

% convert the input XYZ_coords to integer
XYZ_coords=round(XYZ_coords);

% display message
fprintf('searching the area name(s) of the XYZ coordinate(s) for\n');
for ii=1:1:size(XYZ_coords,1)
  fprintf('  %03d: [x,y,z]=[% 4d, % 4d, % 4d]\n',...
          ii,XYZ_coords(ii,1),XYZ_coords(ii,2),XYZ_coords(ii,3));
end
fprintf('in %s space.\n',upper(space));

%% loading VOI atlases
fprintf('setting atlas VOI files to be searched...');
atlas_locations={'automated_anatomical_labeling_VOIs',...
                 'Brainnetome_Atlas_VOIs',...
                 'brodmann_VOIs',...
                 'Buckner_JNeurophysiol11_MNI152_VOIs',...
                 'Choi_JNeurophysiol12_MNI152_VOIs',...
                 'FSL_atlases_VOIs',...
                 'ProbAtlas_v4_VOIs',...
                 'shenparcel_VOIs',...
                 'Talairach.org_VOIs'};

% select the atlases to be used.
if strcmpi(space,'MNI'), atlas_to_be_used(9)=0; end % since 'Talairach.org_VOIs' is defined only in TAL space.
atlas_locations=atlas_locations(logical(atlas_to_be_used));

atlas=cell(length(atlas_locations),1);
if strcmpi(space,'TAL')
  for ii=1:1:length(atlas_locations)
    atlas{ii}=GetFiles(fullfile(fileparts(mfilename('fullpath')),'VOIs',atlas_locations{ii}),'*.voi',{'_TAL';'2mm'});
  end
else % if strcmpi(space,'MNI')
  for ii=1:1:length(atlas_locations)
    atlas{ii}=GetFiles(fullfile(fileparts(mfilename('fullpath')),'VOIs',atlas_locations{ii}),'*.voi',{'_MNI';'2mm'});
  end
end
fprintf('done.\n');

%% searching atlas
roi_name={};
voi_counter=0;
fprintf('searching...\n');
fprintf('\n');
for ii=1:1:length(atlas)
  fprintf('target atlas: %s\n',atlas_locations{ii});
  for jj=1:1:length(atlas{ii})

    % loading the target VOI file
    [voipath,voifname,voiext]=fileparts(atlas{ii}{jj}); %#ok
    fprintf('   searching: %s%s\n',voifname,voiext);
    voi=BVQXfile(atlas{ii}{jj});

    % searching
    find_flg=zeros(size(XYZ_coords,1),1);
    for vv=1:1:voi.NrOfVOIs
      for nn=1:1:size(XYZ_coords,1)
        if ~isempty(intersect(voi.VOI(vv).Voxels,XYZ_coords(nn,:),'rows'))
          find_flg(nn)=1;
          fprintf('              %03d [x,y,z]=[% 4d, % 4d, % 4d] was found in %s.\n',...
                  nn,XYZ_coords(nn,1),XYZ_coords(nn,2),XYZ_coords(nn,3),voi.VOI(vv).Name);
          voi_counter=voi_counter+1;
          roi_name{voi_counter}.name=voi.VOI(vv).Name;          %#ok
          roi_name{voi_counter}.XYZ=XYZ_coords(nn,:);           %#ok
          roi_name{voi_counter}.atlas_type=atlas_locations{ii}; %#ok
          roi_name{voi_counter}.voi_file=atlas{ii}{jj};         %#ok
        end
      end
    end
    for nn=1:1:size(XYZ_coords,1)
      if ~find_flg(nn)
        fprintf('              %03d [x,y,z]=[% 4d, % 4d, % 4d] was not found.\n',...
                nn,XYZ_coords(nn,1),XYZ_coords(nn,2),XYZ_coords(nn,3));
      end
    end

    % clean up
    voi.ClearObject();
    clear voi;
  end % for jj=1:1:length(atlas{ii})
  fprintf('\n');
end % for ii=1:1:length(atlas)
fprintf('completed.\n');

return

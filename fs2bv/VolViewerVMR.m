function VolViewerVMR(VMR_file)

% Loads and displays a BrainVoyager VMR anatomy volume file.
% function VolViewerVMR(VMR_file)
%
% This function loads a BrainVoyager VMR anatomy file and displays it on
% the MATLB figure window with GUI buttons for several functionalities.
%
% [input]
% VMR_file : a BrainVoyager VMR file with a full-path format
%            e.g. VMR_file='F:/Workspace/fMRI_preprocessed/3d/zk11_052.avg.3d_final_TAL.vmr';
%                 or VMR_file='zk11_052.avg.3d_final_TAL.vmr';
%
% [output]
% no output variable, the input VMR_file is displayed on the MATLAB window.
%
% [dependency]
% BVQXtools v0.8d
%
% [reference]
% this function is developed based on imshow3D by Maysam Shahedi,
% distributed on the MATLAB Central File Exchange website.
% https://jp.mathworks.com/matlabcentral/fileexchange/41334-imshow3d--3d-imshow-
%
%
% Created    : "2018-01-26 16:51:39 ban"
% Last Update: "2018-02-22 13:02:01 ban"

% check the input variable
if nargin<1 || isempty(VMR_file), help(mfilename()); return; end

if ~exist(VMR_file,'file'), error('VMR_file not found. check the input variable.'); end

% some constants
global initcoord; % the initial coordinates of the mouse when it is clicked.;
SFntSz=9; % font size

% load the VMR
[vmrpath,vmrfname,vmrext]=fileparts(VMR_file); %#ok
fprintf('loading %s%s...',vmrfname,vmrext);
vmr=BVQXfile(VMR_file);
%vol=rot90(flip(vmr.VMRData,2));
vol=permute(vmr.VMRData,[2,1,3]);
framingcube=vmr.FramingCube;

% adjust the volume size if size(vol) does not match with the framingcube.
xadjust=zeros(3,1);
for ii=1:1:3
  if size(vol,ii)~=framingcube, xadjust(ii)=(framingcube-size(vol,ii))/2/framingcube; end
end

vmr.ClearObject(); clear vmr;
fprintf('done.\n');

% get image intensity
minval=double(min(vol(:)));
maxval=double(max(vol(:)));
cintensity=maxval-minval;
midval=(maxval+minval)/2;
winleveladj=(cintensity+1)/1024;
[Rmin,Rmax]=WL2R(cintensity,midval);

% generate volume display window
figid=figure('Units','pixels','Name',sprintf('BVQX_hbtools VolViewerVMR: %s%s',vmrfname,vmrext),...
             'Color',[0.9,0.9,0.9],'NumberTitle','off','ToolBar','none','MenuBar','none');

% display the input volume with SAG, COR, and TRA views, and set the UI controls
subfig=cell(3,1);
subaxis=cell(3,1);
orienthand=cell(3,1);
stxthand=cell(3,1);
shand=cell(3,1);
orientation_str={'TRA','COR','SAG'};
for vv=1:1:3 % TRA, COR, SAG

  % the subwindow layout below is set by following the BrainVoyager Volume display convention.
  set(figid,'Units','normalized');
  if vv==1 % TRA
    subplot('Position',[0.5+xadjust(3),0.0+xadjust(2),0.5-xadjust(3),0.5-xadjust(2)]); hold on; set(gca,'XLim',[0,size(vol,3)]); set(gca,'YLim',[0,size(vol,2)]);
    subfig{vv}=imshow(squeeze(vol(round(size(vol,vv)/2),:,:)),[minval,maxval]);
  elseif vv==2 % COR
    subplot('Position',[0.5+xadjust(3),0.5+xadjust(1),0.5-xadjust(3),0.5-xadjust(1)]); hold on; set(gca,'XLim',[0,size(vol,3)]); set(gca,'YLim',[0,size(vol,1)]);
    subfig{vv}=imshow(squeeze(vol(:,round(size(vol,vv)/2),:)),[minval,maxval]);
  elseif vv==3 % SAG
    subplot('Position',[0.0+xadjust(1),0.5+xadjust(2),0.5-xadjust(1),0.5-xadjust(2)]); hold on; set(gca,'XLim',[0,size(vol,2)]); set(gca,'YLim',[0,size(vol,1)]);
    subfig{vv}=imshow(squeeze(vol(:,:,round(size(vol,vv)/2))),[minval,maxval]);
  end
  subaxis{vv}=gca;
  set(figid,'Units','pixels');

  % here, set the axis position in pixels, not normalized ratio, just for convenience
  set(subaxis{vv},'Units','pixels');

  % add texts on the image
  axes(subaxis{vv}); %#ok
  orienthand{vv}{1}=text(0.45,0.97,orientation_str{vv},'Color',[1,1,1],'FontWeight','bold','Units','Normalized');
  if vv==1 || vv==2
    orienthand{vv}{2}=text(0.015,0.48,'R','Color',[1,1,1],'FontWeight','bold','Units','Normalized');
    orienthand{vv}{3}=text(0.95,0.48,'L','Color',[1,1,1],'FontWeight','bold','Units','Normalized');
  elseif vv==3
    orienthand{vv}{2}=text(0.015,0.48,'A','Color',[1,1,1],'FontWeight','bold','Units','Normalized');
    orienthand{vv}{3}=text(0.95,0.48,'P','Color',[1,1,1],'FontWeight','bold','Units','Normalized');
  end

  set(get(subaxis{vv},'Children'),'ButtonDownFcn',@mouseClick);

end % for vv=1:1:3 % TRA, COR, SAG

% prepare the positions and sizes of the UI controls
pos_3=get(subaxis{3},'Position');
pos_1=get(subaxis{1},'Position');

base_Pos=[pos_3(1)+35,pos_1(4)-30,80,20];
Labl_Pos=base_Pos+[0,0,0,0];
Hair_Pos=base_Pos+[85,0,15,0];
Wtxt_Pos=base_Pos+[0,-25,-20,0];
Wval_Pos=base_Pos+[65,-25,-40,0];
Ltxt_Pos=base_Pos+[110,-25,-20,0];
Lval_Pos=base_Pos+[110+65,-25,-30,0];

% add UI controls for image intensity controls and decorations
ltxthand=uicontrol('Style','text','Position',Ltxt_Pos,'String','Level','BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz);
wtxthand=uicontrol('Style','text','Position',Wtxt_Pos,'String','Window','BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz);
lvalhand=uicontrol('Style','edit','Position',Lval_Pos,'String',sprintf('%6.0f',midval),'BackgroundColor',[1,1,1],'FontSize',SFntSz,'Callback',@WinLevChanged);
wvalhand=uicontrol('Style','edit','Position',Wval_Pos,'String',sprintf('%6.0f',cintensity),'BackgroundColor',[1,1,1],'FontSize',SFntSz,'Callback',@WinLevChanged);
hairhand=uicontrol('Style','radiobutton','Position',Hair_Pos,'String','Cross-hair','BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz,'Value',1);
lablhand=uicontrol('Style','radiobutton','Position',Labl_Pos,'String','Labels','BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz,'Value',1);

% add sliders
for vv=1:1:3
  S_Pos=get(subaxis{3},'Position'); % set the volume slider under the SAG view
  S_Pos(1)=S_Pos(1)+30; S_Pos(3)=S_Pos(3)-35;
  S_Pos(2)=base_Pos(2)-80; S_Pos(4)=15;

  Stxt_Pos=S_Pos;
  Stxt_Pos(2)=Stxt_Pos(2)+20; Stxt_Pos(4)=15;

  S_Pos=S_Pos+([0,-40,0,0].*(3-vv));
  Stxt_Pos=Stxt_Pos+([0,-40,0,0].*(3-vv));

  % add the UI controls (sliders and slice info strings) and the corresponding callbacks
  if size(vol,vv)>1
    stxthand{vv}=uicontrol('Style','text','Position',Stxt_Pos,...
                          'String',sprintf('%s, Slice# %d / %d',orientation_str{vv},round(size(vol,vv)/2),size(vol,vv)),...
                          'BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz);
  else
    stxthand{vv}=uicontrol('Style','text','Position',Stxt_Pos,...
                          'String','2D image','BackgroundColor',[0.8,0.8,0.8],'FontSize',SFntSz);
  end
  shand{vv}=uicontrol('Style','slider','Min',1,'Max',size(vol,vv),'Value',round(size(vol,vv)/2),...
                      'SliderStep',[1/(size(vol,vv)-1) 10/(size(vol,vv)-1)],'Position',S_Pos,...
                      'Callback',{@SliceSlider,vol,vv});
end

% add the other callback functions
set(figid,'WindowScrollWheelFcn',@mouseScroll);
set(figid,'ButtonDownFcn',@mouseClick);
set(figid,'WindowButtonUpFcn',@mouseRelease);
set(figid,'ResizeFcn',{@figureResized,vol});
set(hairhand,'Callback',@showcrosshair);
set(lablhand,'Callback',@showlabels);

% show cross-hair
hh=[]; hv=[];
showcrosshair(hh,hv);

% here, all positions of the handles are set back to the normalized ratios.
% those translations are important to display images after resizing the MATLAB window.
for vv=1:1:3
  set(subaxis{vv},'Units','normalized');
  set(stxthand{vv},'Units','normalized');
  set(shand{vv},'Units','normalized');
  for ww=1:1:3, set(orienthand{vv}{ww},'Units','normalized'); end
end
set(ltxthand,'Units','normalized');
set(wtxthand,'Units','normalized');
set(lvalhand,'Units','normalized');
set(wvalhand,'Units','normalized');
set(hairhand,'Units','normalized');
set(lablhand,'Units','normalized');

  %% nested functions for the event callbacks

  % figure resize callback
  function figureResized(object,eventdata,vol) %#ok

    % NOTE: here I already tried to control the UI controller sizes in the normalized unit,
    %       however, it resulted in nasty appearance. Therefore, I intentionally use pixel units.
    for pp=1:1:3
      set(subaxis{pp},'Units','pixels');
      set(stxthand{pp},'Units','pixels');
      set(shand{pp},'Units','pixels');
    end
    set(ltxthand,'Units','pixels');
    set(wtxthand,'Units','pixels');
    set(lvalhand,'Units','pixels');
    set(wvalhand,'Units','pixels');
    set(hairhand,'Units','pixels');
    set(lablhand,'Units','pixels');

    pos_3=get(subaxis{3},'Position');
    pos_1=get(subaxis{1},'Position');

    base_Pos=[pos_3(1)+35,pos_1(4)-30,80,20];
    Labl_Pos=base_Pos+[0,0,0,0];
    Hair_Pos=base_Pos+[85,0,15,0];
    Wtxt_Pos=base_Pos+[0,-25,-20,0];
    Wval_Pos=base_Pos+[65,-25,-40,0];
    Ltxt_Pos=base_Pos+[110,-25,-20,0];
    Lval_Pos=base_Pos+[110+65,-25,-30,0];
    set(ltxthand,'Position',Ltxt_Pos);
    set(wtxthand,'Position',Wtxt_Pos);
    set(lvalhand,'Position',Lval_Pos);
    set(wvalhand,'Position',Wval_Pos);
    set(hairhand,'Position',Hair_Pos);
    set(lablhand,'Position',Labl_Pos);

    for pp=1:1:3
      S_Pos=get(subaxis{3},'Position'); % set the volume slider under the SAG view
      S_Pos(1)=S_Pos(1)+30; S_Pos(3)=S_Pos(3)-35;
      S_Pos(2)=base_Pos(2)-80; S_Pos(4)=15;

      Stxt_Pos=S_Pos;
      Stxt_Pos(2)=Stxt_Pos(2)+20; Stxt_Pos(4)=15;

      S_Pos=S_Pos+([0,-40,0,0].*(3-pp));
      Stxt_Pos=Stxt_Pos+([0,-40,0,0].*(3-pp));

      if size(vol,pp)>1
        set(shand{pp},'Position',S_Pos);
      end
      set(stxthand{pp},'Position',Stxt_Pos);
    end

    for pp=1:1:3
      set(subaxis{pp},'Units','normalized');
      set(stxthand{pp},'Units','normalized');
      set(shand{pp},'Units','normalized');
    end
    set(ltxthand,'Units','normalized');
    set(wtxthand,'Units','normalized');
    set(lvalhand,'Units','normalized');
    set(wvalhand,'Units','normalized');
    set(hairhand,'Units','normalized');
    set(lablhand,'Units','normalized');
  end

  % slice slider callback
  function SliceSlider(hObj,event,vol,orientation) %#ok

    S=round(get(hObj,'Value'));
    if orientation==1
      set(subfig{1},'CData',squeeze(vol(S,:,:)));
    elseif orientation==2
      set(subfig{2},'CData',squeeze(vol(:,S,:)));
    elseif orientation==3
      set(subfig{3},'CData',squeeze(vol(:,:,S)));
    else
      error('orientation should be one of 1, 2, and 3. check the input variable.');
    end
    axes(get(subfig{orientation},'Parent'));
    caxis([Rmin,Rmax]);

    if size(vol,orientation)> 1
      set(stxthand{orientation},'String',sprintf('%s, Slice# %d / %d',orientation_str{orientation},S,size(vol,orientation)));
    else
      set(stxthand{orientation},'String','2D image');
    end

    showcrosshair();
  end

  % mouse position in the current figure
  function C=mousePosFig()

    C=get(figid,'CurrentPoint');
  end

  % mouse position in the current axis
  function C=mousePosAxis(axisID)

    C=get(axisID,'CurrentPoint');
  end

  % judge whether the mouse is inside the figure windows
  function [insideX,insideY,orientation]=mouseInsideFig()

    insideX=false; insideY=false;
    orientation=[];
    for pp=1:1:3
      set(subaxis{pp},'Units','pixels');
      P=get(subaxis{pp},'Position');
      C=mousePosFig();
      insideX=C(1,1)<=P(1)+P(3) & C(1,1)>=P(1);
      insideY=C(1,2)>=P(2) & C(1,2)<=P(2)+P(4);
      set(subaxis{pp},'Units','normalized');
      if insideX && insideY
        orientation=pp;
        break
      end
    end
  end

  % judge whether the mouse is inside one of the 3 volume displays
  function [insideX,insideY,orientation]=mouseInsideAxis()

    insideX=false; insideY=false;
    orientation=[];
    for pp=1:1:3
      set(subaxis{pp},'Units','pixels');
      P=get(subaxis{pp},'Position');
      C=mousePosAxis(subaxis{pp});
      C=round([C(1,1),C(1,2)]);
      insideX=0<=C(1) & C(1)<=P(1)+P(3);
      insideY=0<=C(2) & C(2)<=P(2)+P(4);
      set(subaxis{pp},'Units','normalized');
      if insideX && insideY
        orientation=pp;
        break
      end
    end
  end

  % mouse scroll wheel on one of the sub-axes (SAG, COR, TRA image window) callback
  function mouseScroll(object,eventdata) %#ok

    [insideX,insideY,orientation]=mouseInsideFig();
    if insideX && insideY
      UPDN=eventdata.VerticalScrollCount;
      S=get(shand{orientation},'Value')-UPDN;
      S=max(S,1);
      S=min(S,size(vol,orientation));
      if size(vol,orientation)>1
        set(shand{orientation},'Value',S);
        set(stxthand{orientation},'String',sprintf('%s, Slice# %d / %d',orientation_str{orientation},S,size(vol,orientation)));
      else
        set(stxthand{orientation},'String','2D image');
      end

      if orientation==1
        set(subfig{1},'CData',squeeze(vol(S,:,:)));
      elseif orientation==2
        set(subfig{2},'CData',squeeze(vol(:,S,:)));
      elseif orientation==3
        set(subfig{3},'CData',squeeze(vol(:,:,S)));
      else
        error('orientation should be one of 1, 2, and 3. check the input variable.');
      end

      showcrosshair();
    end
  end

  % mouse butoon released callback
  function mouseRelease(object,eventdata) %#ok

    set(gcf,'WindowButtonMotionFcn','');
    %set(gcf,'WindowButtonDownFcn','');
  end

  % mouse click callback
  function mouseClick(object,eventdata) %#ok

    mouseStat=get(gcbf,'SelectionType');

    % window level adjustment
    if mouseStat(1)=='a' % 'alt', mouse right click
      initcoord=get(0,'PointerLocation');
      set(gcf,'WindowButtonDownFcn',@SliceAdj);
      set(gcf,'WindowButtonMotionFcn',@WinLevAdj);

    % cross-hair
    elseif mouseStat(1)=='n' % 'normal', mouse left click
      %initcoord=get(0,'PointerLocation');
      set(gcf,'WindowButtonDownFcn',@SliceAdj);
      set(gcf,'WindowButtonMotionFcn',@SliceAdj);

    end % if mouseStat(1)
  end

  % Window and level mouse adjustment
  function WinLevAdj(varargin)

    if ~isempty(initcoord)
      posdiff=get(0,'PointerLocation')-initcoord;
    else
      posdiff=get(0,'PointerLocation');
    end

    cintensity=cintensity+posdiff(1)*winleveladj;
    midval=midval-posdiff(2)*winleveladj;

    cintensity=min(max(cintensity,1),2*max(vol(:)));
    midval=min(max(midval,1),2*max(vol(:)));

    [Rmin,Rmax]=WL2R(cintensity,midval);
    axes(subaxis{1}); caxis([Rmin,Rmax]);
    axes(subaxis{2}); caxis([Rmin,Rmax]);
    axes(subaxis{3}); caxis([Rmin,Rmax]);
    set(lvalhand,'String',sprintf('%6.0f',midval));
    set(wvalhand,'String',sprintf('%6.0f',cintensity));
    initcoord=get(0,'PointerLocation');

  end

  % window and level text adjustment
  function WinLevChanged(varargin)

    midval=str2double(get(lvalhand,'string'));
    cintensity=str2double(get(wvalhand,'string'));
    if cintensity<1, cintensity=1; end

    [Rmin,Rmax]=WL2R(cintensity,midval);
    axes(subaxis{1}); caxis([Rmin,Rmax]);
    axes(subaxis{2}); caxis([Rmin,Rmax]);
    axes(subaxis{3}); caxis([Rmin,Rmax]);
  end

  % window and level to range conversion
  function [Rmn,Rmx]=WL2R(W,L)

    Rmn=L-(W/2);
    Rmx=L+(W/2);
    if Rmn>=Rmx, Rmx=Rmn+1; end
    %if Rmx>max(vol(:)), Rmx=max(vol(:)); end
  end

  % slice adjusttment by mouse left click and drag
  function SliceAdj(varargin)

    [insideX,insideY,orientation]=mouseInsideAxis();
    if insideX && insideY

      % get the slice value
      pos=mousePosAxis(subaxis{orientation});
      pos=round([pos(1,1),pos(1,2)]);

      if orientation==1 % TRA
        pos(1)=min(max(pos(1),1),size(vol,3));
        pos(2)=min(max(pos(2),1),size(vol,2));

        S(2)=pos(2); % COR
        S(3)=pos(1); % SAG

        set(subfig{2},'CData',squeeze(vol(:,S(2),:)));
        set(subfig{3},'CData',squeeze(vol(:,:,S(3))));
        for pp=[2,3]
          if size(vol,pp)> 1
            set(stxthand{pp},'String',sprintf('%s, Slice# %d / %d',orientation_str{pp},S(pp),size(vol,pp)));
          else
            set(stxthand{pp},'String','2D image');
          end
          set(shand{pp},'Value',S(pp));
        end

      elseif orientation==2 % COR
        pos(1)=min(max(pos(1),1),size(vol,3));
        pos(2)=min(max(pos(2),1),size(vol,1));

        S(1)=pos(2); % TRA
        S(3)=pos(1); % SAG

        set(subfig{1},'CData',squeeze(vol(S(1),:,:)));
        set(subfig{3},'CData',squeeze(vol(:,:,S(3))));
        for pp=[1,3]
          if size(vol,pp)> 1
            set(stxthand{pp},'String',sprintf('%s, Slice# %d / %d',orientation_str{pp},S(pp),size(vol,pp)));
          else
            set(stxthand{pp},'String','2D image');
          end
          set(shand{pp},'Value',S(pp));
        end

      elseif orientation==3 % SAG
        pos(1)=min(max(pos(1),1),size(vol,2));
        pos(2)=min(max(pos(2),1),size(vol,1));

        S(1)=pos(2); % TRA
        S(2)=pos(1); % COR

        set(subfig{1},'CData',squeeze(vol(S(1),:,:)));
        set(subfig{2},'CData',squeeze(vol(:,S(2),:)));
        for pp=[1,2]
          if size(vol,pp)> 1
            set(stxthand{pp},'String',sprintf('%s, Slice# %d / %d',orientation_str{pp},S(pp),size(vol,pp)));
          else
            set(stxthand{pp},'String','2D image');
          end
          set(shand{pp},'Value',S(pp));
        end

      end % if orientation==1 % TRA

      showcrosshair();

    end % if insideX && insideY

  end

  % show cross-hair
  function showcrosshair(object,eventdata) %#ok

    onoff=get(hairhand,'Value');
    if onoff
      S=zeros(3,1);
      for pp=1:1:3, S(pp)=get(shand{pp},'Value'); end

      if isempty(hv) && isempty(hh)
       chcolor=[0,1,1];

        axes(subaxis{1});
        hv(1,1)=plot([S(3),S(3)],[0,S(2)-5],'Color',chcolor,'LineWidth',1.0);
        hv(1,2)=plot([S(3),S(3)],[S(2)+5,size(vol,1)],'Color',chcolor,'LineWidth',1.0);

        hh(1,1)=plot([0,S(3)-5],[S(2),S(2)],'Color',chcolor,'LineWidth',1.0);
        hh(1,2)=plot([S(3)+5,size(vol,3)],[S(2),S(2)],'Color',chcolor,'LineWidth',1.0);

        axes(subaxis{2});
        hv(2,1)=plot([S(3),S(3)],[0,S(1)-5],'Color',chcolor,'LineWidth',1.0);
        hv(2,2)=plot([S(3),S(3)],[S(1)+5,size(vol,1)],'Color',chcolor,'LineWidth',1.0);

        hh(2,1)=plot([0,S(3)-5],[S(1),S(1)],'Color',chcolor,'LineWidth',1.0);
        hh(2,2)=plot([S(3)+5,size(vol,3)],[S(1),S(1)],'Color',chcolor,'LineWidth',1.0);

        axes(subaxis{3});
        hv(3,1)=plot([S(2),S(2)],[0,S(1)-5],'Color',chcolor,'LineWidth',1.0);
        hv(3,2)=plot([S(2),S(2)],[S(1)+5,size(vol,2)],'Color',chcolor,'LineWidth',1.0);

        hh(3,1)=plot([0,S(2)-5],[S(1),S(1)],'Color',chcolor,'LineWidth',1.0);
        hh(3,2)=plot([S(2)+5,size(vol,2)],[S(1),S(1)],'Color',chcolor,'LineWidth',1.0);
      else
        for pp=1:1:3, S(pp)=get(shand{pp},'Value'); end

        set(hv(1,1),'XData',[S(3),S(3)]); set(hv(1,1),'YData',[0,S(2)-5]);
        set(hv(1,2),'XData',[S(3),S(3)]); set(hv(1,2),'YData',[S(2)+5,size(vol,1)]);

        set(hh(1,1),'XData',[0,S(3)-5]); set(hh(1,1),'YData',[S(2),S(2)]);
        set(hh(1,2),'XData',[S(3)+5,size(vol,3)]); set(hh(1,2),'YData',[S(2),S(2)]);

        set(hv(2,1),'XData',[S(3),S(3)]); set(hv(2,1),'YData',[0,S(1)-5]);
        set(hv(2,2),'XData',[S(3),S(3)]); set(hv(2,2),'YData',[S(1)+5,size(vol,1)]);

        set(hh(2,1),'XData',[0,S(3)-5]); set(hh(2,1),'YData',[S(1),S(1)]);
        set(hh(2,2),'XData',[S(3)+5,size(vol,3)]); set(hh(2,2),'YData',[S(1),S(1)]);

        set(hv(3,1),'XData',[S(2),S(2)]); set(hv(3,1),'YData',[0,S(1)-5]);
        set(hv(3,2),'XData',[S(2),S(2)]); set(hv(3,2),'YData',[S(1)+5,size(vol,1)]);

        set(hh(3,1),'XData',[0,S(2)-5]); set(hh(3,1),'YData',[S(1),S(1)]);
        set(hh(3,2),'XData',[S(2)+5,size(vol,2)]); set(hh(3,2),'YData',[S(1),S(1)]);

        drawnow();
      end

      for pp=1:1:3
        for rr=1:1:2
          set(hh(pp,rr),'Visible','on');
          set(hv(pp,rr),'Visible','on');
        end
      end

    else
      for pp=1:1:3
        for rr=1:1:2
          set(hh(pp,rr),'Visible','off');
          set(hv(pp,rr),'Visible','off');
        end
      end
    end % if onoff
  end

  % show orientation and L/R labels
  function showlabels(object,eventdata) %#ok

    onoff=get(lablhand,'Value');
    if onoff
      for pp=1:1:3
        for rr=1:1:3
          set(orienthand{pp}{rr},'Visible','on');
        end
      end
    else
      for pp=1:1:3
        for rr=1:1:3
          set(orienthand{pp}{rr},'Visible','off');
        end
      end
    end
  end

end % function VolViewerVMR(VMR_file)

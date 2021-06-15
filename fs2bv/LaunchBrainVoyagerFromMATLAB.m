function bvqx=LaunchBrainVoyagerFromMATLAB(BV_ver,affinity_mask)

% Launches BrainVoyager QX from MATLAB using Microsoft COM interface.
% function LaunchBrainVoyagerFromMATLAB(BV_ver,:affinity_mask)
% (: is optional)
%
% This function launches BrainVoyager QX from MATLAB
% using Microsoft COM interface
%
% NOTE 1: before using this function and call BrainVoyager from MATLAB,
%         we first need to run SetupBVQX_COM and register BrainVoyager COM
%         object to Windows OS.
% NOTE 2: Windows PowerShell should be installed on your machine.
%
% ******************************************************************************************
% IMPORTANT NOTE : To launch BrainVoyager 20 (ver.3.2) and above, you have to additionally
%                  set the environmental variable, QT_QPA_PLATFORM_PLUGIN_PATH, like
%                  - Variable name : QT_QPA_PLATFORM_PLUGIN_PATH
%                  - Value         : C:\Program Files\BrainVoyager\plugins\platforms
% ******************************************************************************************
%
% [input]
% BV_ver : Major version of BrainVoyager, 1, 2, or 3.
%          **If you are using BV 21 or above, please set this variable to 3**
%          (currently, only 1 and 2 are internally discriminated from the other versions of BV)
%          if empty, this function tries to find the latest BrainVoyager installed
%          and set its version (works if BrainVoyager is installed in the default
%          directory). if the BrainVoyager-installed directory is not found, 3 is set
%          as a default value.
% affinity_mask : (optional) CPU affinity mask that controls which CPUs/threads
%                 to be used, hex value STRING, e.g. affinity_mask='0xFFFFF'
%
% [output]
% bvqx   : COM object to communicate with BrainVoyager
%
%
% Creted     : "2012-07-16 11:02:50 banh"
% Last Update: "2020-03-12 16:50:49 ban"

% check input variable
if nargin<1 || isempty(BV_ver)
  if exist('C:\Program Files\BrainVoyager\BrainVoyager.exe','file')
    BV_ver=3;
  elseif exist('C:\Program Files (x86)\BrainVoyager\BrainVoyagerQX.exe','file')
    BV_ver=2;
  elseif exist('C:\Program Files (x86)\BrainVoyager QX 1.9.9 (32 Bit)\BrainVoyagerQX.exe','file')
    BV_ver=1;
  else
    BV_ver=3;
  end
end

%if nargin<2 || isempty(affinity_mask), affinity_mask=[]; end
if nargin<2 || isempty(affinity_mask)
  machinename=getComputerName();
  if strcmp(machinename,'cvz820')
    if BV_ver<3
      affinity_mask='0xFFFFF';
    else % if BV_ver==3 % or above e.g. BV 20, 21
      affinity_mask=[];
    end
  else
    affinity_mask=[];
  end
end

% initialization
SetupBVQX_COM(BV_ver);

% launch BrainVoyager
try
  if BV_ver==1 % BrainVoyager QX 1.9
    bvqx=actxserver('BrainVoyagerQX.BrainVoyagerQXInterface.1');
  elseif BV_ver==2 % BrainVoyager QX 2.3
    bvqx=actxserver('BrainVoyagerQX.BrainVoyagerQXScriptAccess.1');
  elseif BV_ver==3 % BrainVoyager 20 or above
    if isempty(getenv('QT_QPA_PLATFORM_PLUGIN_PATH'))
      fprintf('setting the environmental variable: QT_QPA_PLATFORM_PLUGIN_PATH to C:\Program Files\BrainVoyager\plugins\platforms...\n');
      setenv('QT_QPA_PLATFORM_PLUGIN_PATH','C:\Program Files\BrainVoyager\plugins\platforms');
    end
    %bvqx=actxserver('BrainVoyagerQX.BrainVoyagerQXScriptAccess.1');
    bvqx=actxserver('BrainVoyager.BrainVoyagerScriptAccess.1'); % from 20.4, the active X server ID has changed
  else
    error('BV_ver should be 1, 2, or 3. check the input variable.')
  end
catch
  fprintf('WARNING: affinity mask can not be set properly. skipping...\n');
end

% set the window size
try % required as BV 20.6 does not support bvqx.ResizeWindow()
  %bvqx.ResizeWindow(800,600);
  bvqx.ResizeWindow(1920,1200);
catch
  % do nothing
end

scrsz=get(0,'ScreenSize');
if scrsz(3)<1920
  sw=round(0.8*scrsz(3));
else
  sw=1920;
end
if scrsz(4)<1200
  sh=round(0.8*scrsz(4));
else
  sh=1200;
end
bvqx.ResizeWindow(sw,sh);
%if BV_ver>2, bvqx.HideLogPane(); end

% set the affinity mask
if ( BV_ver==1 || BV_ver==2 ) && ~isempty(affinity_mask)
  eval(sprintf('!powershell (Get-Process -Name BrainVoyagerQX).ProcessorAffinity=%s',affinity_mask));
elseif BV_ver==3 && ~isempty(affinity_mask)
  eval(sprintf('!powershell (Get-Process -Name BrainVoyager).ProcessorAffinity=%s',affinity_mask));
end

return

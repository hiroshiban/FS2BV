function LaunchBrainVoyagerWithAffinityMask(BV_ver,affinity_mask,BVQX_exe_path)

% Launches BrainVoyager QX from MATLAB as an independent process with a CPU affinity mask
% function LaunchBrainVoyagerWithAffinityMask(BV_ver,:affinity_mask,:BVQX_exe_path)
% (: is optional)
%
% This function launches BrainVoyager QX from MATLAB as an independent process.
% The launched BrainVoyager instance can not be controlled from MATLAB but will
% stay even after the MATLAB process is terminated. An affinity_mask should be
% correctly set with a hex value STRING to control CPUs/threds that will be
% assigned to the BrainVoyager instance.
%
% NOTE 1 : Windows PowerShell should be installed on your machine.
%
% ******************************************************************************************
% IMPORTANT NOTE : To launch BrainVoyager 20 (ver.3.2) and above, you have to additionally
%                  set the environmental variable, QT_QPA_PLATFORM_PLUGIN_PATH, like
%                  - Variable name : QT_QPA_PLATFORM_PLUGIN_PATH
%                  - Value         : C:\Program Files\BrainVoyager\plugins\platforms
% ******************************************************************************************
%
% [input]
% BV_ver        : Major version of BrainVoyager, 1, 2, or 3.
%                 **If you are using BV 21 or above, please set this variable to 3**
%                 (currently, only 1 and 2 are internally discriminated from the other versions of BV)
%                 if empty, this function tries to find the latest BrainVoyager installed
%                 and set its version (works if BrainVoyager is installed in the default
%                 directory). if the BrainVoyager-installed directory is not found, 3 is set
%                 as a default value.
% affinity_mask : (optional) CPU affinity mask that controls which CPUs/threads
%                 to be used, hex value STRING, e.g. affinity_mask='0xFFFFF'
%                 empty by default = use all CPUs/threads
% BVQX_exe_path : (optional) a full path to the BrainVoyagerQX.exe
%                 'C:\Program Files (x86)\BrainVoyager QX 1.9.9 (32 Bit)\BrainVoyagerQX.exe' (when BV_ver==1)
%                 'C:\Program Files (x86)\BrainVoyager\BrainVoyagerQX.exe' (when BV_ver==2)
%                 by default. if you install BrainVoyager somewhere else, please change this value.
%
% [output]
% no output variable,
% the handle etc of the BrainVoyager instance will be displayed on the MATLAB window
%
%
% Creted     : "2015-09-19 14:26:07 ban"
% Last Update: "2018-09-12 12:48:18 ban"

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

if nargin<2 || isempty(affinity_mask), affinity_mask=[]; end
if nargin<3 || isempty(BVQX_exe_path)
  if BV_ver==1
    BVQX_exe_path='C:\Program Files (x86)\BrainVoyager QX 1.9.9 (32 Bit)\BrainVoyagerQX.exe';
  elseif BV_ver==2
    BVQX_exe_path='C:\Program Files (x86)\BrainVoyager\BrainVoyagerQX.exe';
  elseif BV_ver==3 % BV_ver==3 or above % BV 20, 21
    if isempty(getenv('QT_QPA_PLATFORM_PLUGIN_PATH'))
      fprintf('setting the environmental variable: QT_QPA_PLATFORM_PLUGIN_PATH to C:\Program Files\BrainVoyager\plugins\platforms...\n');
      setenv('QT_QPA_PLATFORM_PLUGIN_PATH','C:\Program Files\BrainVoyager\plugins\platforms');
    end
    BVQX_exe_path='C:\Program Files\BrainVoyager\BrainVoyager.exe';
  else
    error('BV_ver should be 1, 2, or 3. check the input variable.')
  end
end

% launch BrainVoyager
try
  if BV_ver==1 || BV_ver==2
    eval(sprintf('!cmd.exe /c start "BrainVoyagerQX.exe" "%s"',BVQX_exe_path));
    if ~isempty(affinity_mask)
      eval(sprintf('!powershell (Get-Process -Name BrainVoyagerQX).ProcessorAffinity=%s',affinity_mask));
    end
  else % if BV_ver==3 % or above e.g. BV 20, 21
    eval(sprintf('!cmd.exe /c start "BrainVoyager.exe" "%s"',BVQX_exe_path));
    if ~isempty(affinity_mask)
      eval(sprintf('!powershell (Get-Process -Name BrainVoyager).ProcessorAffinity=%s',affinity_mask));
    end
  end
catch
  fprintf('WARNING: affinity mask can not be set properly. skipping...\n');
end

return

function [judge,volcorr,volmi,bvvol,bvvolreg,fsvol]=FlipLeftRightCheckerVMR(BV_VMR_file,FS_T1_file,display_flg,thres)

% Check whether the input BrainVoayger VMR anatomy is flipped along the left and right axis.
% function [judge,bvvol,bvvolreg,fsvol]=FlipLeftRightCheckerVMR(BV_VMR_file,FS_T1_file,:display_flg,:thres)
%
% This function loads BrainVoyager VMR anatomy file (*.vmr), compares its orientation
% with a reference FreeSurfer T1 anatomy (eg. T1.mgz), and evaluate whether the VMR is
% unwillingly flipped along the LEFT and RIGHT axis. The VMR's LEFT/RIGHT flip issue is
% always one of the most confusing matters in processing anatomical data obtained without
% any landmark (e.g. a vitamine tablet on the left ear or forehead) unless we can confirm
% the scan orientation or know the details of orientation/position parameters stored in
% the image header. If you are not sure about the image orientations, please check using
% this function.
%
% Internal procedures: to evaluate the BrainVoyager VMR volume flip, this function
% performs the image coregistrations between 1. the 'asis'-VMR (target) and the reference
% FreeSurfer anatomy volume (e.g. T1.mgz or brain.mgz) and 2. 'flipped'-VMR and the
% reference FreeSurfer anatomy. The goodnesses of these two coregistrations are compared
% in terms of correlation-coefficient and mutual information. When the flipped-VMR provides
% a better fit to the reference compared to the asis-VMR, this function judges that the
% taret VMR is flipped along the left and right axis.
%
% [important note]
% Here, the function assumes that the orientation of the reference volume (FreeSurfer T1)
% is correct. please validate it before running this function.
%
% [input]
% BV_VMR_file : BrainVoyager VMR anatomy file (*.vmr), specified with a relative path format,
%               in which the location where this function is called is an origin of the path.
% FS_T1_file  : FreeSurfer anatomy file, T1.mgz, specified with a relative path format,
%               in which the location where this function is called is an origin of the path.
%               This anatomy file is used as a reference and is assumed to be aligned properly
%               along all the superior/inferior, anterior/posterior, and left/right axes.
% display_flg : (optional) whether displaying the coregistration resutl, [0|1].
%               if 1 (non-zero), the coregistration results are displayed. furthermore, the
%               displayed result is saved in the same directory with the input VMR directory
%               as *.png format image. 1 by default.
% thres       : (optional) a threshold used in determining whether the volume is flipped or not.
%               the flip is internally judged by comparing the correlation coeffs or mutual
%               information (MI) between the reference volume (FreeSurfer T1) and VMR (asis vs flip).
%               Then, even when corr_coeff_asis (MI_asis) < corr_coeff_flip (MI_flip) is
%               oberserved, if corr_coeff_flip is lower than the thres, this function returns
%               'unknown', not 'flip' for safety. thres=0.20 by default.
%
% [output]
% judge       : a string, one of 'correct', 'flipped', 'unknown'
%               if the image correlations below are too small, the judge may result
%               in 'unknown'. Then, please manually check the orientation of the data.
% volcorr     : correlation coefficient between 1. asis-VMR vs FreeSurfer_T1, and
%                                               2. flip-VMR vs FreeSurfer_T1
%               if 1 > 2, the input VMR will be in the correct orientation, while
%               if 2 > 1, the input VMR's left and right will be flipped.
% volmi       : mutual information between 1. asis-VMR vs FreeSurfer_T1, and
%                                          2. flip-VMR vs FreeSurfer_T1
%               if 1 > 2, the input VMR will be in the correct orientation, while
%               if 2 > 1, the input VMR's left and right will be flipped.
%               note 1: both volcorr and volmi are taken into account in judging the
%               volume flip. for details, please see the codes below.
%               note 2: volcorr and volmi values would change depending on use_mi_flg variable.
% bvvol       : BrainVoyager raw (uncoregisterd) volume data [x,y,z].
% bvvolreg    : BrainVoyager volume data, coregistered to the reference FreeSurfer anatomy [x,y,z].
% fsvol       : FreeSurfer reference volume data [x,y,z].
%
% [dependency]
% MATLAB Image Processing Toolbox is requried to run this function.
%
%
% Created    : "2017-12-24 10:44:55 ban"
% Last Update: "2018-09-04 18:48:01 ban"

if ~isMATLABToolBoxAvailable('Image Processing Toolbox')
  error('this function requires Image Processing Toolbox. please check the toolbox is available on your machine.');
end

%% check the input variables
if nargin<2, help(mfilename()); return; end
if nargin<3 || isempty(display_flg), display_flg=1; end
if nargin<4 || isempty(thres), thres=0.20; end

if ~exist(fullfile(pwd,BV_VMR_file),'file')
  error('BV_VMR_file not found. check the input variable.');
end

if ~exist(fullfile(pwd,FS_T1_file),'file')
  error('FS_T1_file not found. check the input variable.');
end

%% loading volume datasets

[vmrpath,vmrfname,vmrext]=fileparts(fullfile(pwd,BV_VMR_file));
[fspath,fsfname,fsext]=fileparts(fullfile(pwd,FS_T1_file)); %#ok

% display message
fprintf('checking the target VMR''s Left/Right flip...\n');
fprintf('target (moving)           : %s%s\n',vmrfname,vmrext);
fprintf('reference to be registered: %s%s\n',fsfname,fsext);
fprintf('\n');

% load the BV_VMR_file
fprintf('loading target VRM...');
bvvmr=BVQXfile(fullfile(pwd,BV_VMR_file));
bvvol=bvvmr.VMRData();

% convert the orientation so as to match with the mricron's view when the images are displayed on MATLAB
try
  bvvol=flipdim(permute(bvvol,[2,1,3]),3); %#ok
catch
  bvvol=flip(permute(bvvol,[2,1,3]),3);
end
bvR=imref3d(size(bvvol),bvvmr.VoxResY,bvvmr.VoxResX,bvvmr.VoxResZ);
bvvmr.ClearObject(); clear bvvmr;

fprintf('done.\n');

% load the FS_T1_file
fprintf('converting MGZ to VMR and loading...\n');
ConvertFreeSurferMGZ2VMR(FS_T1_file,'','','','',1,1);
fsvmr=BVQXfile(fullfile(pwd,strrep(FS_T1_file,'.mgz','.mgz.vmr')));
delete(fullfile(pwd,strrep(FS_T1_file,'.mgz','.mgz.vmr'))); % delete the converted VMR file
fsvol=fsvmr.VMRData();

% convert the orientation so as to match with the mricron's view when the images are displayed on MATLAB
try
  fsvol=flipdim(permute(fsvol,[2,1,3]),3); %#ok
catch
  fsvol=flip(permute(fsvol,[2,1,3]),3);
end
fsR=imref3d(size(fsvol),fsvmr.VoxResY,fsvmr.VoxResX,fsvmr.VoxResZ);
fsvmr.ClearObject(); clear fsvmr;

%fprintf('done.\n');

%% prepare coregistration by adjusting the convergence parameters

fprintf('preparing coregistration parameters...');

% [optimizer,metric]=imregconfig('monomodal'); % generate an optimizer and a metrics
optimizer=registration.optimizer.RegularStepGradientDescent;
%metric=registration.metric.MeanSquares; % using mean squared errors
metric=registration.metric.MattesMutualInformation; % using MIs (Mutual Information)

% adjust the optimizer parameter so that the coregistration converges to the global min(max),
% avoiding falling into a local maximum(minimum) for unimodal coregistration
optimizer.GradientMagnitudeTolerance=0.01;
optimizer.MinimumStepLength=0.1;
optimizer.MaximumStepLength=0.5;
optimizer.MaximumIterations=300;
optimizer.RelaxationFactor=0.3;

fprintf('done.\n');

%% run coregistration and judge the image LR flip

fprintf('running coregistration...');
% this is an initial coregistration, if the result is fine, the second one will be skipped.

% run coregistration, bvvol --> fsvox
% here, the 'rigid', not 'affine', transformation should be used for the third parameter
bvvolreg=imregister(bvvol,bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4);

% for validating the image flip, run another coregistration using the flipped VMR
try
  bvvolregflip=imregister(flipdim(bvvol,3),bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4); %#ok
catch
  bvvolregflip=imregister(flip(bvvol,3),bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4);
end

fprintf('done.\n');

% judge whether the input BrainVoayger VMR is flipped or not
fprintf('judging the image flip...');
idx_asis=find(fsvol>0 & bvvolreg>0);
volcorr(1)=corr2(bvvolreg(idx_asis),fsvol(idx_asis)); % correlation between VMR-asis vs reference_FreeSurfer_T1
volmi(1)=mutualInformation(double(bvvolreg(idx_asis)),double(fsvol(idx_asis)));

idx_flip=find(fsvol>0 & bvvolregflip>0);
volcorr(2)=corr2(bvvolregflip(idx_flip),fsvol(idx_flip)); % correlation between VMR-flipped vs reference_FreeSurfer_T1
volmi(2)=mutualInformation(double(bvvolregflip(idx_asis)),double(fsvol(idx_asis)));

fprintf('done.\n');

% run the second coregistration if the initial result is not good
if max(max(volmi),max(volcorr))<thres
  warning off; %#ok
  fprintf('the initial coregistration does not provide a good result.\nrunning the second coregistration...');
  optimizer=registration.optimizer.OnePlusOneEvolutionary;

  % parameer adjustment for multimodal coregistration, e.g. fMRI vs CT
  optimizer.InitialRadius=6.25e-3;%0.001;
  optimizer.Epsilon=1e-3;%1.5e-4;
  optimizer.GrowthFactor=2.0;%2.0;%1.01;
  optimizer.MaximumIterations=300;

  bvvolreg=imregister(bvvol,bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4);

  % for validating the image flip, run another coregistration using the flipped VMR
  try
    bvvolregflip=imregister(flipdim(bvvol,3),bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4); %#ok
  catch
    bvvolregflip=imregister(flip(bvvol,3),bvR,fsvol,fsR,'rigid',optimizer,metric,'PyramidLevels',4);
  end

  fprintf('done.\n');

  % judge whether the input BrainVoayger VMR is flipped or not
  fprintf('judging the image flip...');
  idx_asis=find(fsvol>0 & bvvolreg>0);
  volcorr(1)=corr2(bvvolreg(idx_asis),fsvol(idx_asis)); % correlation between VMR-asis vs reference_FreeSurfer_T1
  volmi(1)=mutualInformation(double(bvvolreg(idx_asis)),double(fsvol(idx_asis)));

  idx_flip=find(fsvol>0 & bvvolregflip>0);
  volcorr(2)=corr2(bvvolregflip(idx_flip),fsvol(idx_flip)); % correlation between VMR-flipped vs reference_FreeSurfer_T1
  volmi(2)=mutualInformation(double(bvvolregflip(idx_asis)),double(fsvol(idx_asis)));

  fprintf('done.\n');
  warning on; %#ok
end

fprintf('\n**********************************************************************\n');
fprintf('Correlations\n');
fprintf('the    ASIS     VMR vs FreeSurfer T1: %.4f\n',volcorr(1));
fprintf('the L/R FLIPPED VMR vs FreeSurfer T1: %.4f\n',volcorr(2));
fprintf('\n');
fprintf('Mutual Information\n');
fprintf('the    ASIS     VMR vs FreeSurfer T1: %.4f\n',volmi(1));
fprintf('the L/R FLIPPED VMR vs FreeSurfer T1: %.4f\n',volmi(2));
fprintf('**********************************************************************\n');

if ( volmi(1)>volmi(2) && volcorr(1)>volcorr(2) ) && ( volmi(1)>thres || volcorr(1)>thres )
  fprintf('the input VMR volume is in the correct orientation.\n');
  judge='correct';
elseif ( volmi(1)<=volmi(2) && volcorr(1)<=volcorr(2) )  && ( volmi(2)>thres || volcorr(2)>thres )
  fprintf('*WARNING* the input VMR volume seems to be FLIPPED. please be careful.\n');
  judge='flipped';
elseif ( volmi(1)>volmi(2) && volcorr(1)<=volcorr(2) ) || ( volmi(1)<=volmi(2) && volcorr(1)>volcorr(2) )
  fprintf('*WARNING* corr-coeff and mutual information are inconsistent. please check by a visual inspection.\n');
  judge='unknown';
else
  fprintf('*WARNING* Too low corr-coeff and mutual information. please check by a visual inspection.\n');
  judge='unknown';
end
fprintf('**********************************************************************\n\n');


%% display the coregistration result

if display_flg
  % set the slice position, following the mricron's convension
  sagpos=129;
  corpos=size(fsvol,2)-115; % as AP is flipped on MATLAB due to the difference of the origin, I have to adjust the slice position as shown here
  trapos=119; %size(fsvol,1)-119; % for legibility, I intentionally changed the slice position for the transverse view.
  scrsz=get(0,'ScreenSize');
  f1=figure('Name',sprintf('BrainVoyager VMR: %s%s to FreeSurfer T1: %s%s. Coregistration result: %s',...
                           vmrfname,vmrext,fsfname,fsext,judge),...
            'Position',[scrsz(3)/7,3*scrsz(4)/8,5*scrsz(3)/7,4*scrsz(4)/8],...
            'Color',[0.9,0.9,0.9],...
            'NumberTitle','off',...
            'ToolBar','none',...
            'Menu','none');

  % sagitgal view
  bigsubplot(1,3,1,1,0.06,0.06);
  imshowpair(fsvol(:,:,sagpos),bvvolreg(:,:,sagpos),'Scaling','joint');
  xlabel('A <-----> P');
  ylabel('I <-----> S');
  title(sprintf('Sag, #slice: %02d',sagpos));
  hold on;
  plot([0,size(fsvol,1)],[trapos,trapos],'-','Color',[1,1,1],'LineWidth',2);
  plot([corpos,corpos],[0,size(fsvol,2)],'-','Color',[1,1,1],'LineWidth',2);

  % coronal view
  bigsubplot(1,3,1,2,0.06,0.06);
  imshowpair(squeeze(fsvol(:,corpos,:)),squeeze(bvvolreg(:,corpos,:)),'Scaling','joint');
  xlabel('L <-----> R');
  ylabel('I <-----> S');
  title(sprintf('Cor, #slice: %02d',size(fsvol,2)-corpos)); % #slice should be adjusted due to the differences of image origin between MATLAB and mricron
  hold on;
  plot([0,size(fsvol,1)],[trapos,trapos],'-','Color',[1,1,1],'LineWidth',2);
  plot([sagpos,sagpos],[0,size(fsvol,3)],'-','Color',[1,1,1],'LineWidth',2);

  % transverse (horizontal) view
  bigsubplot(1,3,1,3,0.06,0.06);
  imshowpair(squeeze(fsvol(trapos,:,:)),squeeze(bvvolreg(trapos,:,:)),'Scaling','joint');
  xlabel('L <-----> R');
  ylabel('P <-----> A');
  title(sprintf('Tra, #slice: %02d',255-trapos)); % #slice should be adjusted due to the differences of image origin between MATLAB and mricron
  hold on;
  plot([sagpos,sagpos],[0,size(fsvol,3)],'-','Color',[1,1,1],'LineWidth',2);
  plot([0,size(fsvol,2)],[corpos,corpos],'-','Color',[1,1,1],'LineWidth',2);

  % saving the image
  set(gcf,'PaperPositionMode','auto');
  print(f1,fullfile(vmrpath,sprintf('%s%s-TO_%s%s_coregistration.png',vmrfname,vmrext,fsfname,fsext)),'-dpng','-r0');
end

return


%% subfunction

% the function below is from mrLoadRet tool.
function [Iab,Pab,Pa,Pb] = mutualInformation(a,b,normalize,nbins)
%
% function [Iab,Pab,Pa,Pb] = mutualInformation(a,b,[nbins])
%
% Computes the mutual information between two vectors. Uses hist2 to
% compute the joint histogram (which ignores the tails of the two marginal
% distributions. Mutual information is:
%    I(a,b) = H(a) + H(b) - H(a,b)
% where
%    H(a) = - sum P(a) log[P(a)]
%
% Normalized mutual information is:
%    [H(a) + H(b)] / H(a,b)
%
% Default nbins: sqrt(length(v1)/10)
% Default normalize: 0
%
% djh, 3/2005

% Default normalize
if ~exist('normalize','var')
    normalize=0;
end

% Default nbins
if ~exist('nbins','var')
    nbins=round(sqrt(length(a)/10));
end

% Joint histogram
abHist=hist2(a,b,nbins);

% Marginal histograms
aHist=sum(abHist,1);
bHist=sum(abHist,2);

% Probabilities
N=sum(aHist);
Pa=aHist/N;
Pb=bHist/N;
Pab=abHist/N;

% Disable divide by 0 and log of 0 warnings
warning('off'); %#ok
Ha=(Pa .* log(Pa));
id=isfinite(Ha);
Ha=- sum(Ha(id));

Hb=(Pb .* log(Pb));
id=isfinite(Hb);
Hb=- sum(Hb(id));

Hab=(Pab .* log(Pab));
id=isfinite(Hab);
Hab=- sum(Hab(id));
warning('on'); %#ok

if normalize
    Iab=(Ha + Hb) / Hab;
else
    Iab=Ha + Hb - Hab;
end

return

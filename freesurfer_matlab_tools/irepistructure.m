function s = irepistructure(TBI,nslices,ndummies,ntp,skip,PreInv,InvDur,TI1,ROFlip)
% s = irepistructure(TBI,nslices,ndummies,ntp,skip,PreInv,InvDur,TI1,ROFlip)

% Pulse sequence parameters
s.TBI = TBI; % 2.092*1000; % time bet inversions in ms, like the TR
s.nslices = nslices;
s.ndummies = ndummies;
s.ntp = ntp;
s.skip = skip; % slice permutation skip, eg 7
s.InvDur = InvDur; % duration of inversion pulse, eg, 12
s.TI1 = TI1; % time from end of inv to first readout
s.ROFlip = ROFlip; % readout flip angle in degrees
s.PreInv = PreInv; % subtract this from TBS of final slice
s.TBS = (s.TBI-s.TI1-s.InvDur+s.PreInv)/s.nslices; % time between slices

% Timing parameters
s.tEvent = [];    % Time of event ms
s.FlipEvent = []; % Degrees, can be 0
s.IsReadOut = []; % 0 or 1 indicating event is a readout event
s.TI = []; % Inversion time at readout
s.EventSliceNo = []; % Anatomical slice number for event
s.AcqSliceNo = [];   % Acquistion slice number for event
s.Slice1PreInv = 0; % Set to 1 to model slice 1 acqed before inv

% Could set timing here with call to irepitiming.m

% MR Parameters
s.T1 = 1200; % msec
s.eff = 1; % efficiency
s.sigma = 0; % noise for rician model
s.biexp = []; % biexponential T1 [fraction T1]

% Cache this so improve speed for rician model
s.lagsnr = [0:.1:10 11:20:10000]';
% This has been commented out because it requires the 
% stats toolbox, and compiled version won't run.
% This data can just be cached and read in too
% This will be used in irepisynth.m
%s.lag  = laguerreL(0.5,-0.5*(s.lagsnr.^2));
if(0)
  % To test
  snr = 100*rand(10,1);
  laghat = interp1(s.lagsnr,s.lag,snr);
  lag = laguerreL(0.5,-0.5*(snr.^2));
end

% Synthesis
s.sliceno = 1; % slice number to synthesize
s.tRead = []; % time of readouts for given slice
s.yMz = []; % magnetization at each Event
s.yRead = []; % magnetization only at readouts

% Optimization
s.parset = 1;

% Fit
s.nexclude = 0; % exclude this number of init slices
s.nminexclude = 0; % exclude this number of mins/nulls
s.iexclude = []; % time points excluded
s.iexclude1 = []; % init time points excluded 
s.iexclude2 = []; % null/min time points excluded 
s.indkeep = []; % indicies of time points included in fit
s.tFit = []; % tRead(s.indkeep)
s.yFit = []; % Data to fit, s.y = s.y0(indkeep,:);
s.y = []; % Data to fit, all time points
s.X0 = []; % Synthesized signal, all time points
s.X = []; % Synthesized signal s.X = s.X0(indkeep,:);, should be XFit
s.M0 = 0; % beta = inv(X'*X)*X'*y;
s.yhat = []; % yhat = X*beta (only includes indkeep)
s.yhat0 = []; % estimate including all time points = beta*x.S0;
s.res = []; % residual y-yhat
s.rstd = 0; % residual variance
s.err = 10e10; % error term, usually std(res)

% Files, etc
s.involfile = '';
s.invol = [];
s.maskfile = '';
s.maskvol = [];
s.configfile = '';
s.outdir = '';

return


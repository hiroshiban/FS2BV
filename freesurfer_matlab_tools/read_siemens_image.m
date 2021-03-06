function m = read_siemens_image(fname)
% m = read_siemens(fname)


%
% read_siemens_image.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_siemens_image: error opening file %s', fname));
	return
end

fseek(fid, 4994, 'bof');
rows = fread(fid, 1, 'short');
cols = fread(fid, 1, 'short');

fseek(fid, 6144, 'bof');
m = fread(fid, [rows cols], 'short');

fclose(fid);

% eof

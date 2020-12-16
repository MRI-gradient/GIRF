function [gOut, kOut, tK] = PredictGrad(this, tIn, gIn, tOut, channels, convType) 
% PredictGrad takes a calculated GIRF and predicts the played out
% trajectory from a nominal input trajectory
%
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


% Assign defaults
filter = 'GIRF'; %{'GIRF' 'PE'}
if ~exist('channels','var')
    channels = {'X' 'Y' 'Z'};
end
if ~exist('convType','var')
    convType = 'conv';
end
if (nargin < 4 ) || isempty(tOut)
    tOut = tIn;
end

% Check dimensions
nCh = length(channels);
nIn = size(gIn,2);
if nIn ~= nCh
    error('Number of inputs (%d) must match number of assigned channels (%s)',nIn,[channels{:}])
end
nOut = size(this.(channels{1}).GIRF, 2);
% nsOut = length(fOut);
nsH = length(this.(channels{1}).f);
nIl = size(gIn,3);

% Assemble H
H = zeros(nsH,nIn,nOut);
for iCh = 1:nCh
    H(:,iCh,:) = this.(channels{iCh}).GIRF;
end
fH = this.(channels{1}).f;
dfH = fH(2)-fH(1);

% Set length of prediction
TH = 1/dfH; % TODO - add option to set TH & dfH for prediction
TO = max(tIn(end),tOut(end)) - min(tIn(1),tOut(1));
TZF = min(TH,TO);

% Prepare input
dtIn = tIn(2) - tIn(1);
nZF = ceil(TZF/dtIn);
gIn = [zeros(nZF,nIn,nIl); gIn; zeros(nZF,nIn,nIl)];
tIn = [tIn(1)+dtIn*(-nZF:-1)'; tIn; tIn(end)+dtIn*(1:nZF)'];
fIn = time2freq(tIn);
IN = fftshift(fft(gIn),1)*dtIn;
nsIn = length(tIn);
T_In = tIn(end) - tIn(1);

% Interpolate GIRF onto input grid
HIp = interp1(fH,H,fIn);
HIp(isnan(HIp)) = 0;

if strcmp(convType,'conv')
    hTime = real(ifft(ifftshift(HIp,1)));
    tStart = 1e-3; % impulse response startup time to use
    tSettle = TZF; % impulse response settling time to use
    tEndInd = ceil(tSettle/dtIn);
    tStartInd = ceil(tStart/dtIn);
    hTime([tEndInd+1:nsIn-tStartInd-1],:,:) = 0;
    HNew = fftshift(fft(hTime),1);
    HNew(isnan(HNew)) = 0;
    HIp = HNew;
end

OUT = zeros(nsIn, nOut, nIl);
for iIl = 1:nIl
    OUT(:,:,iIl) = sum(repmat(IN(:,:,iIl),[1 1 nOut]).*HIp,2);
end
gOut = ifft(ifftshift(OUT,1))/dtIn;
gOut = real(gOut);

% Integrate gradients to k-coefficients
gamma = this.gamma;
kOut = cumsum(gOut)*dtIn*gamma;

% Optionally interpolate output onto new grid
dtOut = tOut(2) - tOut(1);
gOut = interp1(tIn, gOut, tOut);
kOut = interp1(tIn + dtIn/2, kOut, tOut + dtOut/2);
gOut(isnan(gOut)) = 0;
kOut(isnan(kOut)) = 0;

tK = tOut + dtOut/2;






    
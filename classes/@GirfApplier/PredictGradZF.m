function [gOut, kOut, tK] = PredictGradZF(this, tIn, gIn, tOut, predChannels, convType)
%Calculate pre-emphasis prediction
%
% IN
% tIn           Input time vector
% gIn           Input gradients
% tOut          Optional: Input time vector (may be same as input time vector)
% channels      Optional: Channels to predict
% convType      Optional: Type of convolution for prediction
%
% OUT
% gOut          Predicted gradient output (on time vector of tOut)
% kOut          Predicted k-output (on time vector of tOut + dt/2)
% tK            Time vector for predicted k-output
%
% EXAMPLE
%   [gOut, kOut, tK] = girfApplier.PredictGrad(this, tIn, gIn, tOut, channels, convType);
%
%   See also GirfApplier, GirfEssential
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


% Assign default variables 
if ~exist('predChannels','var')
    predChannels = {'X' 'Y' 'Z'};
end
if ~exist('convType','var')
    convType = 'conv';
end
if (nargin < 4 ) || isempty(tOut)
    tOut = tIn;
end

% Check dimensions
nCh = length(predChannels);
nIn = size(gIn,2);
if nIn ~= nCh
    error('Number of inputs (%d) must match number of assigned channels (%s)',nIn,[predChannels{:}])
end
nOut = size(this.girf, 2);
nsH = length(this.freq);
nIl = size(gIn,3);

% Assemble H & select required channels
H = zeros(nsH,nIn,nOut);
for iCh = 1:nCh
    [chExists, chInd] = ismember(predChannels{iCh},this.inChannels);
    if chExists
        H(:,iCh,:) = this.girf(:,:,chInd);
    else
        error('Input channel %s cannot be found',predChannels{iCh})
    end
end
dfH = this.df;
dtH = this.dt;

% Set length of prediction
TH = 1/dfH; 
TO = max(tIn(end),tOut(end)) - min(tIn(1),tOut(1));
% TZF = min(TH,TO);

% Zero-fill input at start and end to avoid aliasing
nsIn = size(gIn,1);
dtIn = tIn(2) - tIn(1);
nZFMin = ceil(TH/dtIn);
[dtR1,dtR2] = rat(dtIn/dtH);
dtN = dtR1*dtR2;
nZFTot = ceil((nZFMin+nsIn)/dtN)*dtN-nsIn;
nZFneg = floor(nZFTot/2);
nZFpos = ceil(nZFTot/2);
gIn = [zeros(nZFneg,nIn,nIl); gIn; zeros(nZFpos,nIn,nIl)];
tIn = [tIn(1)+dtIn*(-nZFneg:-1)'; tIn; tIn(end)+dtIn*(1:nZFpos)'];
fIn = time2freq(tIn);
dfIn = fIn(2)-fIn(1);
IN = fftshift(fft(gIn),1)*dtIn;
nsIn = length(tIn);
T_In = nsIn*dtIn;

% Zero-fill H in time domain to matching T_in
hTime = real(ifft(ifftshift(H,1)));
tStart = 1e-3; % impulse response startup time to use
tStartInd = ceil(tStart/dtIn);
nZFHT = round(T_In/dtH) - size(hTime,1);
hTime = [hTime(1:tStartInd,:,:); zeros(nZFHT,nIn,nOut); hTime(tStartInd+1:end,:,:)];
H = fftshift(fft(hTime),1);
nsH = size(hTime,1);
dfH = 1/(nsH*dtH);
fH = dfH*((1:nsH)' - (floor(nsH/2)+1));

% Zero-fill H or IN to matching bandwidth
if max(fIn) - max(fH) > dfH
    nZFHneg = ceil((fH(1) - fIn(1)))/dfH;
    nZFHpos = ceil((fIn(end)-fH(end)))/dfH;
    H = [zeros(nZFHneg,nIn,nOut); H; zeros(nZFHpos,nIn,nOut)];
    nsH = size(H,1);
    fH = dfH*((1:nsH)' - floor(nsH/2)+1);
    dtH = 1/(nsH*dfH);
elseif max(fH) - max(fIn) > dfH
    dfIn = fIn(2)-fIn(1);
    nZFHneg = ceil((fIn(1) - fH(1)))/dfIn;
    nZFHpos = ceil((fH(end)-fIn(end)))/dfIn;
    IN = [zeros(nZFHneg,nIn,nOut); IN; zeros(nZFHpos,nIn,nOut)];
    nsIn = size(IN,1);
    fIn = dfIn*((1:nsIn)' - floor(nsIn/2)+1);
    dtIn = 1/(nsIn*dfIn);
end
    
% Interpolate H if needed
if abs((dfH-dfIn)/dfH)>1e-6 || abs((dtH-dtIn)/dtH)>1e-6
    H = interp1(fH,H,fIn);
    H(isnan(H)) = 0;
end

OUT = zeros(nsIn, nOut, nIl);
for iIl = 1:nIl
    OUT(:,:,iIl) = sum(repmat(IN(:,:,iIl),[1 1 nOut]).*H,2);
end
gOut = ifft(ifftshift(OUT,1))/dtIn;
gOut = real(gOut);

% Integrate gradients to k-coefficients
gamma = this.gamma1H;
kOut = cumsum(gOut)*dtIn*gamma;

% Optionally interpolate output onto new grid
dtOut = tOut(2) - tOut(1);
gOut = interp1(tIn, gOut, tOut);
kOut = interp1(tIn + dtIn/2, kOut, tOut + dtOut/2);
gOut(isnan(gOut)) = 0;
kOut(isnan(kOut)) = 0;

tK = tOut + dtOut/2;






    
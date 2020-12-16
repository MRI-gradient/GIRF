function PeakElimination(this, fCenter, fWidth, fReps, nFit)
% Interpolates values around artefactual peaks in calculated GIRF
%
% USE
% PeakElimination(fCenter,fWidth,reps)
%
% IN
%   fCenter [Hz] Center of artefactual peak
%   fWidth  [Hz] Width of peak to be interpolated
%   fReps   [Hz] How many peak harmonics to interpolate (default: 1)
%   nFit    Polynomial order of fit (default: 1)
%
% OUT
%   GIRF    [nr_samples nr_k] GIRF after peak interpolation
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


fprintf('Starting peak interpolation...')
ticVS = tic;

% Convert GIRF if not already in frequency domain
if isempty(this.girf) && ~isempty(this.girfTime)
    this.ConvertDomain('time2freq');
end

%% Check input parameters
girf = this.girf;
f = this.freq;
if abs(fCenter) > max(abs(f))
    warning('Peak center frequency too high - no peak interpolation performed')
    return
end
if nargin < 4
    fReps = 1;
end
if nargin < 5
    nFit = 1;
end

%% 
nS = length(f);
nC = floor(nS/2)+1;
for iR = 1:fReps
    indsPeak = find(abs(f-iR*fCenter)<fWidth/2);
    indsPeakNeg = (-length(indsPeak):-1) + 1 + 2*nC-indsPeak(1);
    indsFit = find(abs(f-iR*fCenter)<fWidth*3/2);
    indsFit = setdiff(indsFit,indsPeak);
    nK = size(girf,2);
    for iK = 1:nK
        Preal = polyfit(f(indsFit),real(girf(indsFit,iK)),nFit);
        Pimag = polyfit(f(indsFit),imag(girf(indsFit,iK)),nFit);
        girf(indsPeak,iK) = polyval(Preal,f(indsPeak)) + 1i*polyval(Pimag,f(indsPeak));
    end
    girf(indsPeakNeg,:) = flip(conj(girf(indsPeak,:)));
end

this.girf = girf;
this.ConvertDomain('f2t');

timeVS = toc(ticVS);
fprintf(' done after %f seconds! \n', timeVS)



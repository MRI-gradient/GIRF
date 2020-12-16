function [girf, freq] = ComputeGirf(this)
% Performs the GIRF calculation in the frequency domain
%
%   [girf,freq] = girf.ComputeGirf();
%
% IN
%   timeIn      [nSamples x 1] input time vector 
%   in          [nSamples x nInChannels x nWaveforms] input gradient/shim
%               waveforms
%   timeOut     [nSamples x 1] output time vector (must be same as input)     
%   out         [nSamples x nOutBasis x nWaveforms] measured output
%               gradient/shim waveforms
%
% OUT
%   freq        [nSamples x 1] girf frequency vector 
%   girf        [nSamples x nOutBasis x nInChannels] calculated girf
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


fprintf('Starting GIRF calculation...')
ticGC = tic;

% Check data consistency
if size(this.timeIn,1) ~= size(this.timeOut,1)
    error('Resample input and output waveforms onto same time grid before girf calculation!!')
end

% Get data size
nS          = length(this.timeOut);
nIn         = size(this.in,2);
nOut        = size(this.out,2);
nWaveforms  = size(this.in,3);

% Get frequency domain input and output waveforms
freq    = time2freq(this.timeOut(:,1));
inFreq  = this.inFreq;
outFreq = this.outFreq;

% Find out whether each waveform uses only one inputchannel
if nIn >1
    inWaveforms = zeros(nIn, nWaveforms);
    for inCh = 1: nIn
        for w = 1: nWaveforms
            if numel(find(inFreq(:,inCh,w))) > 4        % non-zero input
                inWaveforms(inCh,w) = 1;                % there is input in this channel and this waveform
            end
        end
    end
useMatrixCalculation = max(sum(inWaveforms,1));         % if more than 1 inputchannel is used in this waveform (max > 1) 
                                                        % matrix calculation needed, else each waveform has only one driving channel
end
% Perform least-squares estimation from inputs
girf = zeros(nS,nOut,nIn);
if nIn == 1
    inSOSInv = 1./sum(abs(inFreq).^2,3);
    for iOut = 1:nOut
        girf(:,iOut) = sum(outFreq(:,iOut,:).*conj(inFreq),3).*inSOSInv;
    end
elseif useMatrixCalculation ==1                         % only one driving channel for each waveform
    for iIn = 1:3
        indexWaveforms  = find(inWaveforms(iIn,:)==1);   % get waveforms which use this input channel
        inSOSInv        = 1./sum(abs(inFreq(:,iIn,indexWaveforms)).^2,3);
        for iOut = 1:nOut
            girf(:,iOut,iIn) = sum(outFreq(:,iOut,indexWaveforms).*conj(inFreq(:,iIn,indexWaveforms)),3).*inSOSInv;
        end
    end   
else                                                    % matrix calculation, more than one inputchannel for one waveform
    for iF = 1:nS
        girf(iF,:,:) = outFreq(iF,:,:)*pinv(inFreq(iF,:,:));
    end
end


girf(isnan(girf)) = 0;

timeGC = toc(ticGC);
fprintf(' done after %f seconds! \n', timeGC)

% Set class properties
this.isFreqDomainGirf   = 1;
this.girf               = girf;
this.freq               = freq;
this.ConvertDomain();                                   % convert to time domain as well
end


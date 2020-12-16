function [girf, freq] = WindowFreq(this, BW, filterType, beta)
% Performs frequency-domain windowing of GIRF
%
% IN
%   girf        [nr_samples nr_k] SIRF matrix to be smoothed
%   freq        [nr_samples 1] frequency vector 
%   BW          [scalar] in [Hz] bandwidth at FWHM of the filter
%   filterType  {'ga' 'bl' 'bh' 'rc'} the filter type (gaussian, blackman, blackman-harris or raised cosine)
%   beta        [scalar] roll-off factor for raised cosine filter
%
% OUT
%   girf        [nr_samples nr_k] smoothed SIRF
%   freq        [nr_samples 1] cut frequency vector
%
% EXAMPLE
%   girf.WindowFreq(this, BW, filterType, beta);
%
%   See also VariableSmoothing, WindowTime
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


% Fill up missing input arguments
if nargin < 3
    filterType = 'rc';
end
if nargin < 4
    beta = 1/3;
end

% Convert GIRF if not already in frequency domain
if isempty(this.girf) && ~isempty(this.girfTime)
    this.ConvertDomain('time2freq');
end

% Perform frequency-domain windowing
[girf, filter] = BW_window(this.girf, this.freq, BW, filterType, beta);
freq = this.freq;

% Cut out frequencies outside of filter
nF = size(filter,1);
inds0 = find(filter(1:floor(nF/2))==0,1,'last');
if ~isempty(inds0)
    girf = girf(inds0:nF+1-inds0,:,:);
    freq = freq(inds0:nF+1-inds0);
end

% Write variables into GIRF object
this.girf = girf;
this.freq = freq;
this.ConvertDomain('freq2time');
this.filterInfo = [this.filterInfo;...
    {sprintf('Freq-domain windowing performed with filter type %s, BW = %d, beta = %d',...
    filterType, BW, beta)}];


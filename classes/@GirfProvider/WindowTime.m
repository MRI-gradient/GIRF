function [girfTime, time] = WindowTime(this, T, filterType, beta)
% Performs time-domain windowing of GIRF
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
%   girf.WindowTime(this, T, filterType, beta);
%
%   See also WindowFreq
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
    beta = 0;
end

% Convert GIRF if not already in time domain
if isempty(this.girfTime) && ~isempty(this.girf)
    this.ConvertDomain('f2t');
end

% Perform time-domain windowing
[girfTime, filter] = BW_window(this.girfTime, this.time, T, filterType, beta);
time = this.time;

% Cut out time samples outside of filter
nS = size(filter,1);
inds0 = find(filter(1:floor(nS/2))==0,1,'last');
if ~isempty(inds0)
    girfTime = girfTime(inds0+1:nS-inds0,:,:);
    time = time(inds0+1:nS-inds0);
end

% Write variables into GIRF object
this.girfTime = girfTime;
this.time = time;
this.ConvertDomain('t2f');
this.filterInfo = [this.filterInfo;...
    {sprintf('Time-domain windowing performed with filter type %s, T = %d, beta = %d',...
    filterType, T, beta)}];


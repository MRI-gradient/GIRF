function [girf, freq] = VarSmoothFreq(this, fMax, BWrange, vsBW)
% Performs variable smoothing on frequency domain GIRF
%
% IN
%   girf        [nr_samples nr_k] SIRF matrix to be smoothed
%   freq        [nr_samples 1] frequency vector 
%   fMax        [Hz] maximum frequency to be smoothed & cut
%   BWrange     [BWmin BWmax] [Hz] minimal and maximal bandwidth of smoothing kernel (optional)
%   vsBW        [Hz] frequency region of narrower smoothing kernel (optional)
%
% OUT
%   girf        [nr_samples nr_k] smoothed SIRF
%   freq        [nr_samples 1] cut frequency vector
%
% EXAMPLE
%   girf.VarSmoothFreq(fMax);
%   girf.VarSmoothFreq(fMax, BWrange, vsBW);
%
%   See also VariableSmoothing
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
    BWrange = [];
end
if nargin < 4
    vsBW = [];
end

% Convert GIRF if not already in frequency domain
if isempty(this.girf) && ~isempty(this.girfTime)
    this.ConvertDomain('time2freq');
end

% Perform variable smoothing
[girf, freq, fMax, BWrange, vsBW] = VariableSmoothing(this.girf, this.freq, fMax, BWrange, vsBW);

% Write variables into GIRF object
this.girf = girf;
this.freq = freq;
this.ConvertDomain('freq2time');
this.filterInfo = [this.filterInfo;...
    {sprintf('Variable smoothing performed with fMax = %d, BWrange = [%d %d], vsBW = %d',...
    fMax, BWrange, vsBW)}];


function ConvertDomain(this,domain)
%Convert frequency-domain GIRF to time domain and vice versa
%
% IN
% domain    Specify conversion domains {'freq2time' 'time2freq'}
%           Default is determined by this.isFreqDomainGirf
%
% OUT
%
% EXAMPLE
%   girf.ConvertDomain();
%   girf.ConvertDomain(freq2time);
%
%   See also GirfEssential
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


if nargin < 2
    if this.isFreqDomainGirf
        domain = 'freq2time';
    else
        domain = 'time2freq';
    end
end
switch domain
    case {'freq2time' 'f2t'}
        % Convert from frequency to time domain
        this.time = time2freq(this.freq);
        this.girfTime = real(fftshift(ifft(ifftshift(this.girf,1)),1));
    case {'time2freq' 't2f'}
        % Convert from time to frequency domain
        this.freq = time2freq(this.time);
        this.girf = fftshift(fft(ifftshift(this.girfTime,1)),1);
end
end

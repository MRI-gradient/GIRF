function [filter] = raised_cosine(f, T, beta)
% creates a raised cosine window for frequency domain filtering
%
% USE
% [filter] = raised_cosine(f, T, beta)
%
% IN
%   f       [nr_samples x 1]; frequency vector 
%   T       [scalar]; 1/BW at FWHM of the window
%   beta    [scalar]; roll-off factor between 0 and 1 for the raised-cosine
%                     window (0 giving a box-car function, and 1 a cosine without plateau)
%
% OUT
%   filter  [nr_samples x 1]; raised cosine window
%
% See also BW_filter
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

filter = zeros(size(f));

ind = abs(f) <= (1+beta)/(2*T);
filter(ind) = (1+cos(pi*T*(abs(f(ind))-(1-beta)/(2*T))/beta))/2;
filter(abs(f) <= (1-beta)/(2*T)) = 1;



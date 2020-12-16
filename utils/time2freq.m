function [f, df, f_max] = time2freq(t)
% Function to compute frequency vector corresponding to a time-vector 
% (assuming fftshift)
%
% USE
% [f, df, f_max] = time2freq(t)
%
% IN
% t         [n_samples x 1] time vector [s]
% 
% OUT
% f         [n_samples x 1] frequency vector [Hz]
% df        frequency resolution [Hz]
% f_max     bandwidth [Hz]
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

nrs = length(t);
dt = t(2)-t(1);
f_max = 1/dt;
df = f_max/nrs;
f = ([0:nrs-1]'-floor(nrs/2))*df;


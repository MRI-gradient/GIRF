function t = CenteredTime(dt, nrs)
% Function to compute time vector, with zero at center (assuming fftshift)
%
% USE
% t = CenteredTime(dt, nrs)
%
% IN
% dt    sampling interval
% nrs   number of samples
% 
% OUT
% t     time vector, centered at zero (assuming fftshift)
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

t = 0:dt:dt*(nrs-1);
cs = floor(nrs/2) + 1;
t = t-t(cs);

function phase = CenteredPhase(array)
% Function to compute unwrapped phase of complex vector, with zero phase at zero
% frequency (assuming fftshift)
%
% USE
% phase = CenteredPhase(array)
%
% IN
% array     [n_samples x n_columns] complex-valued array
% 
% OUT
% phase     [n_samples x n_columns] phase of array
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

phase = unwrap(angle(array));
middle = floor(length(array)/2) + 1;
phase = phase - repmat(phase(middle,:),length(phase),1);

classdef GirfApplier < GirfEssential
%Performs gradient prediction and pre-emphasis computation based on GIRF
%
%
% EXAMPLE
%   girfA = GirfApplier();
%   girfA = GirfApplier(girfEssential);
%
%   See also GirfEssential, GirfProvider
%
% Authors:  Johanna Vannesjo (johanna.vannesjo@gmail.com),
%           Jennifer Nussbaum, Lars Kasper
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


properties
        gamma1H = 2.675222099e8;   % {2.675222099e8} Gyromagnetic ratio of 1H  [rad/s/T]
end % properties
 
 
methods

% Constructor of class
function this = GirfApplier(girfEssential)
    if nargin > 0 && isa(girfEssential,'GirfEssential')
        this.isFreqDomainGirf  = girfEssential.isFreqDomainGirf;
        this.inChannels        = girfEssential.inChannels;
        this.outBasis          = girfEssential.outBasis;

        this.girf              = girfEssential.girf;
        this.freq              = girfEssential.freq;

        this.girfTime          = girfEssential.girfTime;
        this.time              = girfEssential.time;
    end
end

% NOTE: Most of the methods are saved in separate function.m-files in this folder;
%       except: constructor, delete, set/get methods for properties.

end % methods
 
end

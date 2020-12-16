classdef GirfEssential < matlab.mixin.Copyable
%Holds GIRF data in frequency/and or time domain
%
%
% EXAMPLE
%   girf = GirfEssential();
%   girf = GirfEssential(girf,freq,1,inChannels,outBasis);
%   girf = GirfEssential(girfTime,time,0,inChannels,outBasis);
%
%   See also GirfProvider, GirfApplier
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
    % is GIRF calculated in freq or time domain? 
    isFreqDomainGirf; 
    % Input gradient and shim channel names
    inChannels;     
    % Output field basis
    outBasis;       

    % [nSamples x nOutBasis x nInChannels] Calculated GIRF (self- and crossterms)
    girf;           
    % [nSamples] frequency vector for GIRF
    freq;           

    % [nSamples x nOutBasis x nInChannels] GIRF in the time domain
    girfTime;       
    % [nSamples] time vector
    time;           

end % properties
 
properties (Dependent = true, SetAccess = private)
    % [nInChannels] Self-term basis function for each input gradient/shim channel
    selfBasis;      
     
    % nominal frequency resolution of GIRF
    df;            
    % time resolution of GIRF
    dt;            
end % dependent properties
 
methods

% Constructor of class
function this = GirfEssential(data,vect,isFreq,inChannels,outBasis)
    if nargin > 0
        this.isFreqDomainGirf = isFreq;
        if isFreq
            this.girf = data;
            this.freq = vect;
        else
            this.girfTime = data;
            this.time = vect;
        end
        this.ConvertDomain();
        if nargin > 3
            this.inChannels = inChannels;
        end
        if nargin > 4
            this.outBasis = outBasis;
        end
    end
end

% Get methods for dependent properties
function selfBasis = get.selfBasis(this)
    % TODO - use defined enumeration instead!!!
    selfBasis = zeros(1,length(this.inChannels));
    for iIn = 1:length(this.inChannels)
        switch this.inChannels{iIn}
            case {'Z0' 'B0'}
                selfBasis(iIn) = 1;
            case {'X'}
                selfBasis(iIn) = 2;
            case {'Y'}
                selfBasis(iIn) = 3;
            case {'Z'}
                selfBasis(iIn) = 4;
        end
    end
end
function df = get.df(this)
    if isempty(this.freq)
        this.freq = time2freq(this.time);
    end
    df = this.freq(2)-this.freq(1);
end
function dt = get.dt(this)
    if isempty(this.time)
        this.time = time2freq(this.freq);
    end
    dt = this.time(2)-this.time(1);
end

end % methods
 
end

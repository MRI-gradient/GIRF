classdef GirfProvider < GirfEssential
% Class for handling GIRF computation
%
% This class handles all parameteres related to the calculation of a
% measured Gradient or Shim Impulse Response Function.
%
% EXAMPLES: 
% girf = GirfProvider();
% girf = GirfProvider(myGirfEssential);
% girf = GirfProvider(timeIn,in,timeOut,out);
% girf = GirfProvider(timeIn,in,timeOut,out,inChannels,outBasis);
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
        
        %%%% GIRF computation variables %%%%%%%%%%%%%%%%%%%%
        in;             % [nSamples x nInChannels x nWaveforms] Input waveforms to gradient/shim channels in the time domain
        timeIn;         % [nSamples x nWaveforms] time vector for input gradients

        out;            % [nSamples x nOutBasis x nWaveforms] measured gradient/shim fields in the time domain
        timeOut;        % [nSamples x nWaveforms] time vector for measured fields
        
        %%%% Flags and status %%%%%%%%%%%%%%%%
        resampleMethod;

        filterInfo;     % Info about filtering that has been performed

    end
    properties (Dependent = true, SetAccess = private)
        inFreq;         % system inputs in the frequency domain
        outFreq;        % measured gradient/shim fields in the frequency domain

        dtIn;           % gradient input sampling interval
        dtOut;          % sampling interval of measured output      
    end
    methods
        % Class constructor
        function this = GirfProvider(timeIn,in,timeOut,out,inChannels,outBasis)
            if nargin > 0 % Allow for empty constructor
                % Check consistency of data size
                if size(in,3) ~= size(out,3)
                    error('Number of input and output waveforms must match')
                end
                if size(in,1) ~= size(timeIn,1)
                    error('Number of samples in input waveform and input time vector must match')
                end
                if size(out,1) ~= size(timeOut,1)
                    error('Number of samples in output waveform and output time vector must match')
                end

                % Fill up class properties from arguments
                this.timeIn = timeIn;
                this.in = in;
                this.timeOut = timeOut;
                this.out = out;
                
                % Optionally fill up channel information properties
                if nargin > 4
                    if length(inChannels) ~= size(in,2)
                        warning('Different number of input channels specified in input waveform and input channel information')
                    end
                    this.inChannels = inChannels;
                end
                if nargin > 5
                    if length(outBasis) ~= size(out,2)
                        warning('Different number of output basis terms specified in output waveform and output basis information')
                    end
                    this.outBasis = outBasis;
                end
%                 if nargin > 6
%                     if length(selfBasis) ~= size(in,2) || max(selfBasis) > size(out,2)
%                         warning('Check consistency of self-term basis info!')
%                     end
%                     this.selfBasis = selfBasis;
%                 end
            end
        end % Constructor
        
        % Get methods for dependent properties
        function inFreq = get.inFreq(this)
            inFreq = fftshift(fft(this.in),1);
        end
        function outFreq = get.outFreq(this)
            outFreq = fftshift(fft(this.out),1);
        end
        function dtIn = get.dtIn(this)
            dtIn = this.timeIn(2)-this.timeIn(1);
        end
        function dtOut = get.dtOut(this)
            dtOut = this.timeOut(2)-this.timeOut(1);
        end
        
        % Return GirfEssential object
        function girfEssential = GetGirfEssential(this)
            girfEssential                   = GirfEssential();
            
            girfEssential.isFreqDomainGirf  = this.isFreqDomainGirf;
            girfEssential.inChannels        = this.inChannels;
            girfEssential.outBasis          = this.outBasis;
            
            girfEssential.girf              = this.girf;
            girfEssential.freq              = this.freq;
            
            girfEssential.girfTime          = this.girfTime;
            girfEssential.time              = this.time;
        end
    end
end
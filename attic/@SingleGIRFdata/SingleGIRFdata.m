classdef SingleGIRFdata < matlab.mixin.Copyable
% Class for handling GIRFs and trajectory predictions
%
% girfX = SingleGIRFdata();
% 
% This class handles all parameteres related to the calculation of a
% measured Gradient or Shim Impulse Response Function. 
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%
    
    
    properties
        
        %%%% GIRF - essential variables %%%%%%%%%%%%%%%%%%%%
        GIRF;           % Calculated GIRF (self- and crossterms)
        f;              % frequency vector for GIRF
        df;             % nominal frequency resolution of GIRF

        inputs;         % [nSamples x nIn] System inputs in the time domain
        outputs;        % [nSamples x nIn x nOut] measured gradient/shim fields in the time domain
        nIn;            % Number of input pulses
        nOut;           % Number of output field basis functions
        
        tIn;            % time vector for input gradients
        dtIn;           % gradient input sampling interval
        tOut;           % time vector for measured fields
        dtOut;          % sampling interval of measured output      
        
        %%%% Switches and filter parameters %%%%%%%%%%%%%%%%
        doZerofill = [];% [s] Fill up in/out to 1 Hz frequency resolution
        doVarSmoothing = 0; % Perform variable smoothing on calculated GIRF
        isSmoothed = 0; % Has variable smoothing been performed already?

        windowType = 'rc';% Type of BW window {'ga' 'rc' 'bl' 'bh'}
        BW = [];           % BW of GIRF window, FWHM
        BW_rc = 1/3;      % Transition parameter of BW window (only for raised cosine window)
        isWindowed = 0; % Has BW windowing been performed already?

        self;           % Gradient or shim channel monitored field number
        channel;        % Gradient or shim channel name
        
        %%%% Input properties (all) %%%%%%%%%%%%%%%%%%%%%
        pulseType       % input pulse type {'blips' 'sweeps' 'arb'}
        ipType = 'ideal';    %{'ideal' 'sampleAndHold' 'interpolated'}
                
        %%%% Blip input related properties %%%%%%%%%%%%%%
        dur;            % input blip slope durations
        dur2;           % input blip slope durations
        plateau;        % input trapezoid plateau lengths
        fixedSlope = 0; % pulse slope if same for all blips
        amp;            % input pulse amplitudes
        t0;             % input pulse starting time (relative to tOut)

        %%%% Sweep input related properties %%%%%%%%%%%%%%
        sweepType;      % sweep type {'linear' 'quadratic' 'cubic' int}
        Tsweep;         % sweep pulse length
        ampSweep;       % input pulse amplitudes
        t0Sweep;        % input pulse starting time (relative to tOut)
        f1;             % sweep start frequency
        f2;             % sweep end frequency
        phi0;           % sweep pulse starting phase
        AM;             % structure determining sweep amplitude modulation
        smoothing;      % structure deteriming sweep end smoothing
        slew = 200;     % slew rate limit of gradient system
        ilPhaseShift = 0;% shift of sweep pulse between interleaves
        ilTshift = 0;    % shift factor of Tsweep between interleaves

        %%%% Measured output data handling - obsolete?? %%%%%%%%%%%%%%%
%         gamma = 2.675222099e8; % 1H, as in ReconstructionData

        dataPath;       % path to reconstructed traj-/or NI-data
        origDataPath = 'data';  % path to original dataset (for finding params folder)
        dataID;         % [n_concats x nIn] data numbers to be loaded
        kBasis;         % selected output field basis functions  TODO - is this used at all???
        dyns;           % {nIn} dynamics to be averaged
        il = 1;         % interleaves to use for GIRF calculation
        ilTot = 1;      % total number of interleaves acquired (for correct shifting)
        tUse;           % time interval of measured data to use for GIRF calculation
        
        doOffsetCorr;  % correct for baseline offset {0 or 1}
        doSlopeCorr;   % correct for baseline slope {0 or 1}
        tCorr              % time interval for calculation of correction parameters
                        
        %%%% Obsolete or useless...? %%%%%%%%%%%%%%%%%%%%%%%%
        isReadinInput = 0;   % If input acquired with Philips readin patch
        nr_averages;
        average_k;
        TR;
        shifts;
        Grads;
                
    end
    properties (Dependent = true, SetAccess = private)
        IN;     % system inputs in the frequency domain
        OUT;    % measured gradient/shim fields in the frequency domain
        t_k;    % time vector shifted by dt/2 (for phase data)
        girft;  % GIRF in the time domain
    end
    methods

        % Get methods for dependent properties
        function IN = get.IN(this)
            IN = fftshift(fft(this.inputs),1);
        end
        function OUT = get.OUT(this)
            OUT = fftshift(fft(this.outputs),1);
        end
        function t_k = get.t_k(this)
            t_k = this.tOut+this.dt/2;
        end
        function girft = get.girft(this)
            girft = fftshift(ifft(ifftshift(this.GIRF,1)),1);
        end
        
    end
end
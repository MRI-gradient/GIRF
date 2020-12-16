classdef SystemGIRFdata < matlab.mixin.Copyable & dynamicprops
% Class for handling GIRFs and trajectory predictions
%
% girfo = SystemGIRFdata();
% 
% This class contains calculated GIRFs/SIRFs for all channels of a system, 
% and handles GIRF-based output prediction, pre-emphasis calculation
% and pre-emphasis filtering. 
%
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

    
    properties (Dependent = true, SetAccess = private)
        GIRF;   
        f;      
        channels; % channels with calculated GIRFs
    end
    
    properties
        
        %%%% Different single GIRFs %%%%%%%%%%%%%%%%%%%%%%%%%
        Z0
        X
        Y
        Z
        XY
        ZY
        Z2
        ZX
        X2_Y2
        Y3
        XYZ
        Z2Y
        Z3
        Z2X
        ZX2_ZY2
        X3
        Z4
        
        allChannels = {'Z0','X','Y','Z','XY','ZY','Z2','ZX','X2_Y2', ...
                'Y3','XYZ','Z2Y','Z3','Z2X','ZX2_ZY2','X3','Z4'};
                        
        %%% Preemphasis related parameters %%%%%%%%%%%%%%%%%%
        PE;     % Calculated preemphasis filter
        PEType = 'self'; % Set PE caclulation flags {'self' 'matrix'}
        PEInChannels = {}; % Abstract input channels to use
        PESysChannels = {}; % Physical input channels to use
        PEOutChannels = {}; % Output channels to use
        PESetDC; % Flag for how to handle deviant DC responses {'normalize' 'set1' ''}
        Ht;     % Target system response
        HtType = {'erf'}; % Type of target system response filter {'erf' 'rc'}
        HtBW = 30e3;   % Bandwidth of target filter, FWHM
        HtBeta; % Roll-off factor for raised cosine target filter

        %%% Basic physical constants %%%%%%%%%%%%%%%%%%%%%%%%
        gamma = 2.675222099e8; % 1H, as in ReconstructionData
    end

    methods      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function this = SystemGIRFdata(varargin)

            % Initiate SingleGIRFdata objects
            this.Z0     = SingleGIRFdata();
            this.X      = SingleGIRFdata();
            this.Y      = SingleGIRFdata();
            this.Z      = SingleGIRFdata();
            this.XY     = SingleGIRFdata();
            this.ZY     = SingleGIRFdata();
            this.Z2     = SingleGIRFdata();
            this.ZX     = SingleGIRFdata();
            this.X2_Y2  = SingleGIRFdata();
            this.Y3     = SingleGIRFdata();
            this.XYZ    = SingleGIRFdata();
            this.Z2Y    = SingleGIRFdata();
            this.Z3     = SingleGIRFdata();
            this.Z2X    = SingleGIRFdata();
            this.ZX2_ZY2 = SingleGIRFdata();
            this.X3     = SingleGIRFdata();
            this.Z4     = SingleGIRFdata();
            
            % Set self-terms of SingleGIRFdata objects
            for iG = 1:16
                this.(this.allChannels{iG}).self = iG;
                this.(this.allChannels{iG}).channel = this.allChannels{iG};
            end
            this.Z4.self = 7;
            this.Z4.channel = 'Z4';
        end % constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Return dependent variables %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Return calculated GIRFs from SingleGIRFdata objects
        function girf = get.GIRF(this)
            ch = this.channels;
            for iG = 1:length(ch)
                girf(:,:,iG) = this.(ch{iG}).GIRF;
            end
        end
        
        % Return frequency vector of calculated GIRFs
        function f = get.f(this)
            ch = this.channels;
            f = this.(ch{1}).f;
        end

        % Return channels with calculated GIRFs
        function ch = get.channels(this)
            j = 1;
            for iG = 1:17
                if ~isempty(this.(this.allChannels{iG}).GIRF)
                    ch{j} = this.allChannels{iG};
                    j = j+1;
                end
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Save, Load, Clear functions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Save GIRF data either as object, or just GIRF and f
        function Save(this,saveName,varargin)
            if nargin > 2
                saveType = varargin{1};
            else
                saveType = 'variables';
            end
            switch saveType
                case {0 'obj' 'object'}
                    girfo = this;
                    save(saveName, 'girfo')
                case {1 'var' 'variables'}
                    GIRF = this.GIRF;
                    f = this.f;
                    channels = this.channels;
                    save(saveName, 'GIRF', 'f', 'channels')
                case {2 'dir' 'directories'}
                    if ~exist(saveName,'dir')
                        mkdir(saveName)
                    end
                    channels = this.channels;
                    for iCh = 1:length(channels)
                        if ~isempty(this.GIRF{iCh})
                            fileName = [saveName filesep channels{iCh}];
                            GIRF = this.GIRF{iCh};
                            f = this.f{iCh};
                            channel = channels{iCh};
                            save(fileName, 'GIRF', 'f', 'channel')
                        end
                    end
            end
        end

        % Load GIRF data saved as variables
        function Load(this,loadName,loadType,varargin)
            if nargin < 3
                if exist(loadName,'var')
                    loadType = 'var';
                elseif exist(loadName,'dir')
                    loadType = 'dir';
                else
                    error('Could not recognize %s as neither file nor directory',loadName)
                end
            end
            
            switch loadType
                case {'var','variables'}
                    try
                        GIRFdata = load(loadName);
                    catch
                        error('Could not load: %s', loadName)
                    end
                    
                    if isfield(GIRFdata,'channels')
                        channels = GIRFdata.channels;
                    else
                        channels = this.channels;
                    end
                    for iCh = 1:length(channels)
                        this.(channels{iCh}).GIRF = GIRFdata.GIRF(:,:,iCh);
                        this.(channels{iCh}).f = GIRFdata.f;
                    end
                case {'dir','directories'}
                    for iCh = 1:length(this.allChannels)
                        ch = this.allChannels{iCh};
                        fileName = [loadName filesep ch '.mat'];
                        if exist(fileName,'file')
                            try
                                GIRFdata = load(fileName);
                                if isfield(GIRFdata,'SIRF')
                                    this.(ch).GIRF = GIRFdata.SIRF;
                                else
                                    this.(ch).GIRF = GIRFdata.GIRF;
                                end
                                this.(ch).f = GIRFdata.f;
                            catch
                                fprintf('Could not find SIRF/GIRF variable in: %s \n',fileName)
                            end
                        else
                            fprintf('Could not find file: %s \n',fileName)
                        end
                    end
            end
        end
        
        % Remove inputs & outputs from single GIRF data to save disk space
        function ClearInOut(this)
            for iG = 1:17
                this.(this.channels{iG}).inputs = [];
                this.(this.channels{iG}).outputs = [];
            end            
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

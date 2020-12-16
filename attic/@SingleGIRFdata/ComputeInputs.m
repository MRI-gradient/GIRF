function inputs = ComputeInputs(this)
% Function to compute set of blip inputs
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


% Check type of input pulse
if isempty(this.pulseType) || ~iscell(this.pulseType)
    ex1 = '    girfo.pulseType = {''blips''}';
    ex2 = '    girfo.pulseType = {''sweeps''}';
    ex3 = '    girfo.pulseType = {''blips'' ''sweeps''}';
    error('Please specify input pulse type as cell object, examples: \n %s \n %s \n %s',ex1,ex2,ex3)
end

inputsAll = [];
for iP = 1:length(this.pulseType)
    switch this.pulseType{iP};
        case {'blips' 'blip' 'Blips' 'Blip'}
            isBlips = 1; isSweeps = 0;
        case {'sweeps' 'sweep' 'Sweeps' 'Sweep'}
            isSweeps = 1; isBlips = 0;
    end
    
    % Set blip parameters
    if isBlips
        nPulses = length(this.dur);
        if this.fixedSlope
            this.amp = this.fixedSlope*this.dur;
        end
        if length(this.t0)<nPulses
            this.t0 = this.t0*ones(size(this.dur));
        end
        if length(this.amp)<nPulses
            this.amp = this.amp*ones(size(this.dur));
        end
        if isempty(this.plateau)
            this.plateau = zeros(size(this.dur));
        end
        if isempty(this.dur2)
            this.dur2 = this.dur;
        end
    end
    
    % Set sweep parameters
    if isSweeps
        if iscell(this.Tsweep)
            nPulses = length(this.Tsweep);
        else
            nPulses = 1;
            this.Tsweep = {this.Tsweep};
            this.ampSweep = {this.ampSweep};
            this.phi0 = {this.phi0};
            this.f1 = {this.f1};
            this.f2 = {this.f2};
            this.sweepType = {this.sweepType};
            this.t0Sweep = {this.t0Sweep};
            this.slew = {this.slew};
            this.AM = {this.AM};
            this.smoothing = {this.smoothing};
        end
    end
    
    % Compute input pulse
    if strcmp(this.ipType,'ideal')
        if isBlips
            inputs = trapezoid(this.tOut, this.t0, this.amp, this.dur, this.plateau, this.dur2);
        elseif isSweeps
            for iIn = 1:nPulses
                inputs(:,iIn) = sweeps(this.tOut, this.Tsweep{iIn}, this.f1{iIn}, ...
                    this.f2{iIn}, this.phi0{iIn}, this.ampSweep{iIn}, this.t0Sweep{iIn}, ...
                    this.sweepType{iIn}, this.AM{iIn}, this.smoothing{iIn}, this.slew{iIn});
            end
        end
    else
        inputs = zeros(length(this.tOut),nPulses);
        for iIn = 1:nPulses
            if isBlips
                t0 = this.t0(iIn);
                pulseDur = this.dur(iIn)+this.plateau(iIn)+this.dur2(iIn);
                tIn = (0:this.dtIn:pulseDur+this.dtIn)-this.dtIn/2;
                pulseIn = trapezoid(tIn, 0, this.amp(iIn), this.dur(iIn), this.plateau(iIn), this.dur2(iIn));
            elseif isSweeps
                t0 = this.t0Sweep{iIn};
                pulseDur = sum(this.Tsweep{iP});
                tIn = (0:this.dtIn:pulseDur+this.dtIn)-this.dtIn/2;
                pulseIn = sweeps(tIn, this.Tsweep{iIn}, this.f1{iIn}, ...
                    this.f2{iIn}, this.phi0{iIn}, this.ampSweep{iIn}, 0, ...
                    this.sweepType{iIn}, this.AM{iIn}, this.smoothing{iIn}, this.slew{iIn});
            end
            
            if strcmp(this.ipType,'sampleAndHold') % discretize on input raster
                indsOut = find(this.tOut >= t0 & this.tOut <= t0+pulseDur);
                indsIn = ceil((this.tOut(indsOut)-t0)/this.dtIn)+1;
                inputs(indsOut,iIn) = pulseIn(indsIn);
            elseif strcmp(this.ipType,'interpolated') % calculate on input raster & do linear interpolation onto output raster
                indsOut = find(this.tOut >= tIn(1)+ t0 & this.tOut <= tIn(end)+ t0);
                inputs(indsOut,iIn) = interp1(tIn+t0,pulseIn,this.tOut(indsOut));
            end
        end
    end
    % Concatenate all input pulses
    inputsAll = [inputsAll inputs];
end

this.inputs = inputsAll;
this.nIn = size(inputsAll,2);

end



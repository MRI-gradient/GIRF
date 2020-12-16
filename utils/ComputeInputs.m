function inputs = ComputeInputs(dtIn, tOut, params, pulseType, ipType)
% Function to compute set of blip inputs
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
if isempty(pulseType) || ~iscell(pulseType)
    ex1 = '    girfo.pulseType = {''blips''}';
    ex2 = '    girfo.pulseType = {''sweeps''}';
    ex3 = '    girfo.pulseType = {''blips'' ''sweeps''}';
    error('Please specify input pulse type as cell object, examples: \n %s \n %s \n %s',ex1,ex2,ex3)
end

inputsAll = [];
for iP = 1:length(pulseType)
    switch pulseType{iP};
        case {'blips' 'blip' 'Blips' 'Blip'}
            isBlips = 1; isSweeps = 0;
        case {'sweeps' 'sweep' 'Sweeps' 'Sweep'}
            isSweeps = 1; isBlips = 0;
    end
    
    % Set blip parameters
    if isBlips
        nPulses = length(params.dur);
        if params.fixedSlope
            params.amp = params.fixedSlope*params.dur;
        end
        if length(params.t0)<nPulses
            params.t0 = params.t0*ones(size(params.dur));
        end
        if length(params.amp)<nPulses
            params.amp = params.amp*ones(size(params.dur));
        end
        if isempty(params.plateau)
            params.plateau = zeros(size(params.dur));
        end
        if isempty(params.dur2)
            params.dur2 = params.dur;
        end
    end
    
    % Set sweep parameters
    if isSweeps
        if iscell(params.Tsweep)
            nPulses = length(params.Tsweep);
        else
            nPulses = 1;
            params.Tsweep = {params.Tsweep};
            params.ampSweep = {params.ampSweep};
            params.phi0 = {params.phi0};
            params.f1 = {params.f1};
            params.f2 = {params.f2};
            params.sweepType = {params.sweepType};
            params.t0Sweep = {params.t0Sweep};
            params.slew = {params.slew};
            params.AM = {params.AM};
            params.smoothing = {params.smoothing};
        end
    end
    
    % Compute input pulse
    if strcmp(ipType,'ideal')
        if isBlips
            inputs = trapezoid(tOut, params.t0, params.amp, params.dur, params.plateau, params.dur2);
        elseif isSweeps
            for iIn = 1:nPulses
                inputs(:,iIn) = sweeps(tOut, params.Tsweep{iIn}, params.f1{iIn}, ...
                    params.f2{iIn}, params.phi0{iIn}, params.ampSweep{iIn}, params.t0Sweep{iIn}, ...
                    params.sweepType{iIn}, params.AM{iIn}, params.smoothing{iIn}, params.slew{iIn});
            end
        end
    else
        inputs = zeros(length(tOut),nPulses);
        for iIn = 1:nPulses
            if isBlips
                t0 = params.t0(iIn);
                pulseDur = params.dur(iIn)+params.plateau(iIn)+params.dur2(iIn);
                tIn = (0:dtIn:pulseDur+dtIn)-dtIn/2;
                pulseIn = trapezoid(tIn, 0, params.amp(iIn), params.dur(iIn), params.plateau(iIn), params.dur2(iIn));
            elseif isSweeps
                t0 = params.t0Sweep{iIn};
                pulseDur = sum(params.Tsweep{iP});
                tIn = (0:dtIn:pulseDur+dtIn)-dtIn/2;
                pulseIn = sweeps(tIn, params.Tsweep{iIn}, params.f1{iIn}, ...
                    params.f2{iIn}, params.phi0{iIn}, params.ampSweep{iIn}, 0, ...
                    params.sweepType{iIn}, params.AM{iIn}, params.smoothing{iIn}, params.slew{iIn});
            end
            
            if strcmp(ipType,'sampleAndHold') % discretize on input raster
                indsOut = find(tOut >= t0 & tOut <= t0+pulseDur);
                indsIn = ceil((tOut(indsOut)-t0)/dtIn)+1;
                inputs(indsOut,iIn) = pulseIn(indsIn);
            elseif strcmp(ipType,'interpolated') % calculate on input raster & do linear interpolation onto output raster
                indsOut = find(tOut >= tIn(1)+ t0 & tOut <= tIn(end)+ t0);
                inputs(indsOut,iIn) = interp1(tIn+t0,pulseIn,tOut(indsOut));
            end
        end
    end
    % Concatenate all input pulses
    inputsAll = [inputsAll inputs];
end

inputs = inputsAll;

end



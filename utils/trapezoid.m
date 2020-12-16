function [grads] = trapezoid(t,ons,amp,dur,varargin)
% Function to create trapezoidal pulse shapes from parameters
%
% USE
% grads = trapezoid(t,ons,amp,dur,plateau,dur2)
%
% IN
% t         [n_samples x 1] time vector [s]
% ons       [1 x n_pulses] onset of pulse [s] 
% amp       [1 x n_pulses] amplitude of pulse
% dur       [1 x n_pulses] duration of first slope
% plateau   [1 x n_puless] duration of pulse plateau (opt: varargin{1})
% dur2      [1 x n_pulses] duration of second slope (opt: varargin{2})
% 
% OUT
% grads     [n_samples x n_pulses] trapezoidal gradient pulse
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


% Initialize parameters
n_pulses = length(amp); % number of pulses
ns       = length(t);   % number of samples
if size(t,1)<size(t,2), t = t'; end
grads = zeros(ns, n_pulses);
if nargin < 6
    dur2 = dur;  % both slopes equal, if not otherwise specified
else
    dur2 = varargin{2};
end
if nargin < 5
    plateau = zeros(1,n_pulses); % no plateau, if not otherwise specified
else
    plateau = varargin{1};
end

% Create gradient pulses
for iP = 1:n_pulses
    
    % Find timing indices
    t0 = find(t>ons(iP),1,'first');
    ts1 = find(t<ons(iP)+dur(iP),1,'last');
    if plateau(iP) ~= 0
        tp = find(t<ons(iP)+dur(iP)+plateau(iP),1,'last');
    end
    ts2 = find(t<ons(iP)+dur(iP)+plateau(iP)+dur2(iP),1,'last');
    
    slope1 = amp(iP)/dur(iP);
    slope2 = amp(iP)/dur2(iP);
    
    grads(t0:ts1,iP) = slope1*(t(t0:ts1)-ons(iP));
    grads(ts1+1:ts2,iP) = amp(iP) - slope2*(t(ts1+1:ts2)-(ons(iP)+plateau(iP)+dur(iP)));
    if plateau(iP) ~= 0
        grads(ts1+1:tp,iP) = amp(iP);
    end
    
end


function [sweep, f_t, phiEnd] = sweeps(t, T_acq, f1, f2, phi0, A, t_start, type, varargin)
% Function to create frequency-swept pulses from parameters
%
% USE
% sweep = sweeps(t, T_acq, f1, f2, phi0, A, t_start, type)
%
% IN
% t         [n_samples x 1] time vector [s]
% T_acq     length of sweep [s] 
% f1        starting frequency
% f2        end frequency
% phi0      starting phase
% A         pulse amplitude
% t_start   pulse start time
% type      which frequency-modulation of pulse?
% 
% OUT
% sweep     [n_samples x 1] frequency-swept pulse
% f_t       [n_samples x 1] frequency-modulation as function of time
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

% assign optional input variables
if nargin > 8
    AM = varargin{1};
else
    AM.type = 'none';
    AM.slew = [];
end
if nargin > 9
    smoothing = varargin{2};
else
    smoothing.type = 'none';
end
if nargin > 10
    slew = varargin{3};
else
    slew = [];
end

% Initialize sweep pulse
% sweep = sin(phi0)*ones(size(t));
sweep = zeros(length(t),1);
f_t = zeros(length(t),1);

% Check if there is a wish to concatenate sweeps
if iscell(type)
    n_segments = length(type);
    if n_segments > 1
        [sweep1, f_t1, phiEnd1] = sweeps(t, T_acq(1:end-1), f1, f2(1:end-1), phi0, 1, t_start, {type{(1:end-1)}});
        sweep = sweep + sweep1;
        f_t = f_t + f_t1;
        f1 = f2(end-1);
        f2 = f2(end);
        phi0 = phiEnd1;
        t_start = t_start + sum(T_acq(1:end-1));
        T_acq = T_acq(end);
        type = type{end};
    else
        type = type{1};
    end
end

% Single sweep calculation
if ~ischar(type)
    n = type;
else
    % Set parameters according to old chirp script
    switch type
        case {'beta' 'linear'}
            n = 2;
        case 'gamma'
            n = 3;
        case 'delta'
            n = 4;
        case 'betaAM'
            n = 2;
            AM.type = 'gaussian';
            AM.A_max = 10;
            AM.width = 1;
        case 'gammaAM'
            n = 3;
            AM.type = 'gaussian';
            AM.A_max = 10;
            AM.width = 1;
        case 'deltaAM'
            n = 4;
            AM.type = 'gaussian';
            AM.A_max = 10;
            AM.width = 1;
    end
end

% shift time vector and find indices inside sweep
t = t-t_start;
t_inds = find(t>=0 & t<=T_acq);

% calculate sweep
beta = 2*pi*(f2-f1)./(n*T_acq.^(n-1));
sweep(t_inds) = sin( beta*t(t_inds).^n + 2*pi*f1*t(t_inds) + phi0);
f_t(t_inds) = (n*beta/(2*pi))*t(t_inds).^(n-1) + f1;
phiEnd = mod(beta*T_acq^n + 2*pi*f1*T_acq + phi0, 2*pi);
% end

% Scale sweep with possibly frequency-dependent amplitude
AM_boost = ones(size(sweep));
switch AM.type
    case 'gaussian'
        AM_boost(t_inds) = (AM.A_max-1)*exp(-((f_t(t_inds)/AM.width).^2)) + 1; 
end

AM_smooth = ones(size(sweep));
switch smoothing.type
    case 'erf'
        cf = f2 - smoothing.width/2;
        AM_smooth(t_inds) = 0.5 - 0.5*erf((f_t(t_inds)-cf)/(smoothing.width/3)); % Amplitude modulation for smoothing
    case 'time-defined'
        smooth_time = smoothing.width;
        smooth_ct = T_acq - smooth_time/2;
        AM_smooth(t_inds) = 0.5 - 0.5*erf((t(t_inds)-smooth_ct)/(smooth_time/3)); % Amplitude modulation for smoothing
end
 
A_t = A*AM_boost;
if slew
    S_limit = abs(slew./(2*pi*f_t));
    S_limit(isnan(S_limit)) = 0;
    A_t(t_inds) = sign(A)*min(abs(A_t(t_inds)),S_limit(t_inds));
end

sweep = A_t.*AM_smooth.*sweep;
% sweep = A_t.*sweep;
if size(sweep,1) < size(sweep,2)
    sweep = sweep';
end






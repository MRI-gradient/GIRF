% Script to test GIRF code
%
% Authors:  Johanna Vannesjo (johanna.vannesjo@gmail.com),
%           Lars Kasper
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

clear all, clear classes

%% Compute input pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create input time vector
dtIn = 10e-6;
tIn = (0:dtIn:500e-3)';

% Define input sweep
sweepType = 'linear';
Tsweep = 50e-3;
f1 = 0;
f2 = 40e3;
phi0 = 0;
A = 1;
tStart = 5e-3;
inSweep = sweeps(tIn,Tsweep,f1,f2,phi0,A,tStart,sweepType);

% Define input blips
A = 10;
slope = 0.1e-3;
plateau = 0;
tStart = 5e-3;
inTrapz = trapezoid(tIn,tStart,A,slope,plateau);

in(:,1,1) = inSweep;
in(:,1,2) = inTrapz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulate output pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
dtOut = dtIn;
tOut = tIn;
outBasis = 1:16;
out = zeros(length(tOut),length(outBasis),size(in,3));

% Simulate system low-pass filtering
selfOut = BW_filter(in,tIn,20e3,'rc',0.6);
crossOut = 0.1*BW_filter(in,tIn,10e3,'rc',1);
out(:,4,:) = selfOut;
out(:,1,:) = crossOut;

% Add simulated noise to integrated waveform
outK = dtOut*cumsum(out);
outK = outK + 1e-7*randn(size(outK));
out = [zeros(1,length(outBasis),size(in,3)); diff(outK)/dtOut];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute GIRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialize GIRF object %%%%%%
inChannels = {'Z'};
girfo = GirfProvider(tIn, in, tOut, out, inChannels, outBasis);
girfo.Vis('inout','t');
girfo.Vis('inout','f');
girfo.Vis('sensitivity');

%%%% Compute GIRF %%%%%%%%%%%%%%%%%%%%%%%
girfo.ComputeGirf();
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');
girfo.Vis('GIRF','f',[],1);

%%%% Filter GIRF %%%%%%%%%%%%%%%%%%
girfo.WindowFreq(60e3,'rc',1/4);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

girfo.VarSmoothFreq(30e3);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

girfo.WindowTime(10e-3,'rc',0);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

%%%% Return GirfEssential object %%%%
girfE = girfo.GetGirfEssential();

%%%% Test perform a prediction %%%%
girfA = GirfApplier(girfE);
[gOut, kOut, tK] = girfA.PredictGrad(tIn, in(:,:,1), tOut, inChannels, 'conv');
figure, plot(tIn, in(:,:,1), tOut, out(:,4,1), tOut, gOut(:,4))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('girf was calculated')


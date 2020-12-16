
function [k_data] = BW_filter(k_data, t, BW, filterType, beta)
% BW_filter applies a frequency domain filter to the data
%
% The Function BW_filter applies a frequency-domain filter to the data.
% Typical input-data are k-evolutions that can be filtered given the finite
% bandwith of the gradient chain.
% The filter BW corresponds to the FWHM. Three different filters are
% disposable. The raised cosine filter is set as default. However, the filter BW makes much more % of a difference than
% the filter-type.
% 
% USE
% [k_data] = BW_filter(k_data, t, BW, filterType, beta)
%
% IN
%   k_data      [nr_samples nr_orders nDim3]
%   t           [nr_samples 1] in [sec]
%   BW          [scalar] in [Hz] bandwidth at FWHM of the filter
%   filterType  [char], optional: 'ga', 'bl', 'bh' or 'rc': the filter type (gaussian, blackman, blackman-harris or raised cosine)
%   beta_rc     [scalar], optional: roll-off factor between 0 and 1 for the raised cosine filter (0 giving a box-car function, and 1 a cosine without plateau)
%
% OUT
%   k_data      [nr_samples x nr_orders]; filtered data
%
% Authors:   Johanna Vannesjo (johanna.vannesjo@gmail.com), 
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%



%% Return if no filtering desired
if isempty(BW) || BW==0
    return
end

%% internal input
% use the raised-cosine filter if not otherwise specified
if nargin == 3
    filterType = 'rc';
    beta = 1/3;
end 

%% check and adjust format of input data
if ismatrix(k_data) && size(k_data,2)>size(k_data,1)
    k_data = k_data.';
end
if size(t,2)>size(t,1)
    t = t.';
end
if size(t,1) ~= size(k_data,1)
    error('Number of samples in the input data and time vector must agree')
end

if max(max(abs(imag(k_data))))==0
    is_real = 1;
else
    is_real = 0;
end

%% shift data to fourier domain center by subtractig its linear phase
%% (haeb)
if ~is_real;
    linear_phase = zeros(size(k_data));
    for cnt1 = 1:size(k_data,2)
        for cnt2 = 1:size(k_data,3)
            tmp = squeeze(unwrap(angle(k_data(:,cnt1,cnt2))));
            b = regress(tmp, [ones(size(t(:,1))) t(:,1)]);
%             omega(:,cnt1,cnt2) = b(2);
            linear_phase(:,cnt1,cnt2) = b(2).*t(:,1)';
        end
    end
    k_data = k_data .* exp(-1i*linear_phase);
end

%% make data an even function for DFT (haeb)
k_data = [flipud(k_data); k_data];

%% transform data to fourier domain
k_data = ifftshift(fft(fftshift(k_data,1)),1);
dt = t(2)-t(1);
nrs = size(k_data,1);
df = 1/(nrs*dt);
f = ([0:nrs-1]'-floor(nrs/2))*df;

%% apply the filter in frequency space
k_data = BW_window(k_data, f, BW, filterType, beta);

%% transform data back to time domain
k_data = ifftshift(ifft(fftshift(k_data,1)),1);

%% remove symmetric part again
k_data = k_data(size(k_data,1)/2+1:end,:,:);

%% undo phase removal
if ~is_real, k_data = k_data .* exp(1i*linear_phase); end

%% remove falsely introduced imaginary parts
if is_real, k_data = real(k_data); end




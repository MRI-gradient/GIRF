function [array, filter] = BW_window(array, f, BW, filterType, beta)
% Applies selected bandwidth window to input array
%
% The Function BW_window applies a frequency-domain windowing to the data.
% The filter BW corresponds to the FWHM. Available filters are Gaussian,
% Blackman, Blackmanharris and Raised cosine. 
%
% IN:
%   array       [nr_samples x nr_orders] data to be windowed
%   f           [nr_samples x 1] in [Hz] frequency vector of the array
%   BW          [scalar] in [Hz] bandwidth at FWHM of the filter
%   filterType  {'ga' 'bl' 'bh' 'rc'} the filter type (gaussian, blackman, blackman-harris or raised cosine)
%   beta        [scalar] roll-off factor for raised cosine filter
%
% OUT:
%   array       [nr_samples x nr_orders] windowed data
%   filter      [nr_samples x 1] calculated filter
%
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

%% Internal input
% use the Blackman-Harris filter if not otherwise specified
if nargin == 3
    filterType = 'bh';
end 

%% check and adjust format of input data
if ismatrix(array) && size(array,2)>size(array,1)
    array = array.';
end

if size(f,2)>size(f,1)
    f = f.';
end
if size(f,1) ~= size(array,1)
    error('Number of samples in the input data and frequency vector must agree')
end

%% compute the frequency filter
df = f(2)-f(1);
cs = floor(length(f)/2)+1;
filter = zeros(size(f));
switch filterType
    case 'ga' %%% Gaussian filter
        filter = exp(-((f/(BW/2)).^2)*log(2));
        %disp(['-> a gaussian filter is used'])
    case 'rc' %%% Raised cosine filter
        filter = raised_cosine(f,1/BW,beta);
    case 'bl' %%% Blackman filter [FWHM = 405.1 out of 1000]
        ns = floor(BW/df*1000/405.1/2); 
        b = blackman(2*ns+1); % force odd number of samples
        if length(b) <= length(f)
            filter(cs-ns:cs+ns) = b;
        else
            bi = (1:length(f))-cs+ns+1;
            filter = b(bi);
            warning('the chosen filter-BW is too wide for the Blackman filter')
        end
        %disp(['-> a Blackman filter is used'])
    case 'bh' %%% Blackman-Harris filter [FWHM = 342.8 out of 1000]
        ns = floor(BW/df*1000/342.8/2);
        bh = blackmanharris(2*ns+1);  % force odd number of samples
        if length(bh) <= length(f)
            filter(cs-ns:cs+ns) = bh;
        else
            bi = (1:length(f))-cs+ns+1;
            filter = bh(bi);
            warning('the chosen filter-BW is too wide for the Blackman-Harris filter')
        end
        %disp(['-> a Blackman-Harris filter is used'])
end

%% apply the filter in frequency space
array = array.*repmat(filter,[1 size(array,2) size(array,3)]);



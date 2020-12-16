function [SIRF, fs, fMax, BWrange, vsBW] = VariableSmoothing(SIRF, f, fMax, BWrange, vsBW)
% Performs frequency-dependent smoothing of GIRF/SIRF
%
% USE
% [SIRF] = VariableSmoothing(SIRF, f, f2)
%
% IN
%   SIRF    [nr_samples nr_k] SIRF matrix to be smoothed
%   f       [nr_samples 1] frequency vector 
%   fMax    [Hz] maximum frequency to be smoothed & cut
%   BWrange [BWmin BWmax] [Hz] minimal and maximal bandwidth of smoothing kernel (optional)
%   vsBW    [Hz] frequency region of narrower smoothing kernel (optional)
%
% OUT
%   SIRF    [nr_samples nr_k] smoothed SIRF
%   fs      [nr_samples 1] cut frequency vector
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

fprintf('Starting variable smoothing...')
ticVS = tic;

%% Check input parameters
if ~exist('BWrange','var') || isempty(BWrange)
    BW_Kmin = 0.05; % minimal bandwidth of smoothing kernel
    BW_Kmax = 10; % maximal bandwidth of smoothing kernel
else
    BW_Kmin = BWrange(1);
    BW_Kmax = BWrange(2);
end
if ~exist('vsBW','var') || isempty(vsBW)
    vsBW = 40; % frequency region of narrower smoothing kernel
end

%% Create smoothing kernels (K)

% Create time and frequency vector of kernel
df = f(2)-f(1);
Ts = 10/df; % sets oversampled frequency resolution of K (10 * T)
dtk = 0.01/BW_Kmax; % sets oversampled time resolution of k (0.01 * dt)
tk = CenteredTime(dtk,Ts/dtk); %-pb/2;
fK = time2freq(tk);

% Create parameters for time-domain window
pb = 1/BW_Kmax;   % flattop length in seconds
beta = 1/3;       % width of transition band
TK = pb/(1-beta); % length of kernel
BW = 1/TK;        % FWHM bandwidth of kernel

% Create kernel of widest bandwidth
k = raised_cosine(tk, BW, beta).';
K = fftshift(fft(ifftshift(k)));
K = K./max(abs(K)); % normalize in frequency domain

% Cut freq-domain side lobes of kernel
nSL = 3; % number of side lobes to include
ind0 = find(abs(K)<1e-6); % find samples close to 0
cs = ceil(length(ind0)/2); 
K = K(ind0(cs-nSL):ind0(cs+nSL+1));
fK = fK(ind0(cs-nSL):ind0(cs+nSL+1));

% Select SIRF frequency region to smooth
csS = floor(size(SIRF,1)/2)+1; % center sample of SIRF
nrf1 = fix(-fK(1)/df);  % half-width of K within SIRF frequency vector
nrf2 = fix(fK(end)/df); % half-width of K within SIRF frequency vector
nrBW = floor(fMax/df);  % SIRF region to smooth
ind_SK = csS-nrf1:csS+nrBW+nrf2; % indices to smooth (only positive frequencies)
f = f(ind_SK);
SIRF = SIRF(ind_SK,:,:);

% Create temporary to-be-smoothed SIRF vector
nS = nrBW+1;
nOut = size(SIRF, 2);
nIn = size(SIRF,3);
S_temp = zeros(nS,nOut,nIn);

% Create scaled versions of kernel, to perform variable smoothing
ft = f(1:1+nrf1+nrf2);
vsr = BW_Kmax/BW_Kmin; % variable smoothing ratio
sci = find(f<vsBW, 1, 'last'); % SIRF center indices to be smoothed with freq-dependent kernel
scf = (1/vsr-1)*exp(-(f(1+nrf1:sci)/20).^2)+1; % Gaussian scaling of freq-dependent bandwidth
% scf = (1-1/vsr)*f(1+nrf1:sci)/f(sci) + 1/vsr; % linear scaling
% scf = -(vsr-1)*(1-cos(pi*f(1+nrf1:sci)/vsBW))/2 + vsr; % cosine scaling
nrscf = length(scf);    % number of separate smoothing kernels
Kmat = zeros(length(ft),nS);
for l = 1:nrscf
    Ki = interp1(fK*scf(l),K,ft); % interpolate scaled version of smoothing kernel
    Ki(isnan(Ki)) = 0;
    Ki = Ki.*exp(1i*2*pi*ft/scf(l)*pb./2); % shift kernel to center at start of passband
    Ki = Ki./abs(sum(Ki));      % normalize interpolated smoothing kernel
    Kmat(:,l) = Ki;             % keep all kernels in one matrix
end
Ki = interp1(fK,K,ft);          % create unscaled kernel
Ki = Ki.*exp(1i*2*pi*ft*pb/2);  % shift kernel to center at start of passband
Ki = Ki./abs(sum(Ki));
Kmat(:,nrscf+1:end) = repmat(Ki,1,nS-nrscf);

% Smooth SIRF with frequency-dependent kernel, K
for j = 1:nS
    S_temp(j,:,:) = sum(SIRF(j:j+nrf1+nrf2,:,:).*repmat(Kmat(:,j),[1 nOut nIn]));
end

% Mirror smoothed SIRF to negative frequencies
SIRF = [conj(S_temp(end:-1:2,:,:)); S_temp];
fs = [-f(nrf1+1+nrBW:-1:nrf1+2); f(nrf1+1:nrf1+1+nrBW)];

timeVS = toc(ticVS);
fprintf(' done after %f seconds! \n', timeVS)



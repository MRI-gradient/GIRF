function [ PE ] = ComputePE( this )
% Compute preemphasis filters from measured GIRFs
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


fprintf('Starting PE filter calculation...')
ticPE = tic;

%% Load & scale data

% Determine channels for which to calculate PE filter
if ~isempty(this.PESysChannels) % Physical system channels
    sysCh = this.PESysChannels;
else
    sysCh = this.channels;
    this.PESysChannels = sysCh;
end
nSys = length(sysCh);

if ~isempty(this.PEInChannels) % Abstract input channels
    inCh = this.PEInChannels;
else
    inCh = sysCh;
    this.PEInChannels = inCh;
end
nIn = length(inCh);

if ~isempty(this.PEOutChannels) % Output channels
    outCh = this.PEOutChannels;
    nOut = length(outCh);
    for iO = 1:nOut
        out(iO) = this.(outCh{iO}).self;
    end
else
    nOut = size(this.(sysCh{1}).GIRF, 2);
    out = 1:nOut;
end

% Assemble H
f = this.(sysCh{1}).f;
df = f(2)-f(1);
nsH = length(f);
H = zeros(nOut,nSys,nsH);
DC_scaling = 1;
for iCh = 1:nSys
    ch = sysCh{iCh};
    if strcmp(this.PESetDC,'normalize')
        DC_scaling = 1/this.(ch).GIRF(f==0,this.(ch).self);
    end
%     SIRF_all(:,ind,:) = DC_scaling*diag(Hz_scaling)*(SIRF(:,shim_select).');
    H(:,iCh,:) = DC_scaling*this.(ch).GIRF(:,out).';
    if strcmp(this.PESetDC,'set1')
        H(out==this.(ch).self,iCh,f==0) = 1;
    end
end

%% Define desired transfer function, Ht

nHt = length(this.HtBW);
if nHt ~= 1 && nHt ~= nIn
    error('The number of selected target bandwidths must either be 1 or match the number of abstract input channels')
end
Ht = zeros(length(f),nHt);
for iH = 1:nHt
    switch this.HtType{iH}
        case 'erf'
            % Define step response from error function
            dt = 1/(nsH*df);
            t = CenteredTime(dt,nsH);
            srf = (erf(2*this.HtBW(iH)*t)+1)/2;
            irf = [0 diff(srf)];
            Ht(:,iH) = abs(fftshift(fft(ifftshift(irf))))';
        case 'rc'
            % Use raised cosine functions
            Ht(:,iH) = raised_cosine(f, 1/this.HtBW(iH), this.HtBeta(iH));
    end
end
this.Ht = Ht;
if nHt == 1
    Ht = repmat(Ht,[1 nSys]);
end

%% Perform matrix inversion to determine pre-emphasis

PE = zeros(nSys, nIn, nsH);
if strcmp(this.PEType,'self')
    for iS = 1:nSys
        PE(iS,iS,:) = Ht(:,iS)./permute(H(out==this.(sysCh{iS}).self,iS,:),[3 1 2]);
    end
else
    HtTemp = zeros(nOut,nIn,nsH);
    for iI = 1:nIn
        HtTemp(out==this.(inCh{iI}).self,iI,:) = Ht(:,iI).';
    end
    Ht = HtTemp;
    for iF = 1:nsH
        %     Hp_all(:,:,indf) = pinv(SIRF_all(:,:,indf))*W(indf)*diag(Hz_ref(shim_select));
        PE(:,:,iF) = pinv(H(:,:,iF))*Ht(:,:,iF);
    end
end
this.PE = PE;

timePE = toc(ticPE);
fprintf(' done after %f seconds! \n', timePE)

end


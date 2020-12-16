function GIRF = ComputeGirf(this)
% Actual GIRF computation
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


fprintf('Starting GIRF calculation...')
ticGC = tic;

% Get data size
ns = length(this.tOut);
if isempty(this.nOut)
    this.nOut = size(this.outputs,3);
end
if isempty(this.nIn)
    this.nIn = size(this.inputs,2);
end

% Call utils resample function to get input and output onto same raster!
% Need flag for using input or output time raster as base for calculations?
% 

% Zero-fill input and output for increased GIRF frequency resolution
if ~isempty(this.doZerofill)
    dt = this.tOut(2)-this.tOut(1);
    T = ns*dt;
    ns_fill = ceil((this.doZerofill - T)/dt);
    inputs = [this.inputs; zeros(ns_fill, this.nIn)];
    IN = fftshift(fft(inputs),1);
    outputs = zeros(ns+ns_fill, this.nIn, this.nOut);
    outputs(1:ns,:,:) = this.outputs;
    OUT = fftshift(fft(outputs),1);
    tOut = [0:dt:(ns+ns_fill-1)*dt]' + this.tOut(1);
    ns = ns+ns_fill;
else
    IN = this.IN;
    OUT = this.OUT;
    tOut = this.tOut;
end

% Perform least-squares estimation from inputs
f = time2freq(tOut);
IN_inv = 1./sum(abs(IN).^2,2);
GIRF = zeros(ns,this.nOut);
for iK = 1:this.nOut
    GIRF(:,iK) = IN_inv.*sum(conj(IN).*OUT(:,:,iK),2);
end
GIRF(isnan(GIRF)) = 0;

timeGC = toc(ticGC);
fprintf(' done after %f seconds! \n', timeGC)

% if this.doVarSmoothing
%     [GIRF, f] = VariableSmoothing(GIRF, f, this.BW);
% else
%     GIRF = BWwindow(GIRF, f, this.BW, this.windowType, this.BW_rc);
% end

this.GIRF = GIRF;
this.f = f;
end


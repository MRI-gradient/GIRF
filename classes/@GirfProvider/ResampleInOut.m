function ResampleInOut( this, timeRef )
%Resamples input and output gradients onto reference time waveform
%
%   output = ResampleInOut(timeRef)
%
% IN
%
% OUT
%
% EXAMPLE
%   ResampleInOut
%
%   See also
%
% Author:   Jennifer Nussbaum (nussbaum@biomed.ee.ethz.ch)
% Copyright (C) 2017 IBT, University of Zurich and ETH Zurich
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


%% resample Input
if size(timeRef,1)~=size(this.in,1)
    x       = zeros(size(timeRef,1),size(this.in,2),size(this.in,3));
    
    for n1 = 1:size(this.in(1:size(this.in,1),:),2)
        x(:,n1) = interp1(this.timeIn,this.in(:,n1),timeRef);
    end
    
    this.timeIn=timeRef;
    this.in=x;
end
%% resample Output
if size(timeRef,1)~=size(this.out,1)
    x       = zeros(size(timeRef,1),size(this.out,2),size(this.out,3));
    
    for n1 = 1:size(this.out(1:size(this.out,1),:),2)
        x(:,n1) = interp1(this.timeOut,this.out(:,n1),timeRef);
        
    end
    
    this.timeOut=timeRef;
    this.out=x;
end
            
end
function [ GIRF, f ] = FilterGIRF( this )
%Window & smooth GIRf
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


if this.doVarSmoothing
    if ~this.isSmoothed
        [GIRF, f] = VariableSmoothing(this.GIRF, this.f, this.BW/2);
        this.isSmoothed = 1;
        this.f = f;
    end
else
    if ~this.isWindowed
        GIRF = BWwindow(this.GIRF, this.f, this.BW, this.windowType, this.BW_rc);
        this.isWindowed = 1;
    end
end
this.GIRF = GIRF;

end


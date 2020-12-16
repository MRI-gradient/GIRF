function inputs = load_input_grads(this)
% Function to load input gradients from files in parameter
% directory
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


this.nIn = size(this.dataID,2);

% Load from file
for iN = 1:this.nIn
    [G, t_in] = load_grads_from_file(1, [], 1, this.dataID(1,iN));
    %                 in_temp(:,iN) = G(:,this.self-1);
    in_temp(:,iN) = G(:,2);
end

% Interpolate onto output time grid
ns_in = length(t_in);
ns_out = length(this.tOut);
dt_in = t_in(2) - t_in(1);
dtOut = this.tOut(2) - this.tOut(1);
if ns_in ~= ns_out || dt_in ~= dtOut
    for iN = 1:this.nIn
        inputs(:,iN) = interp1(t_in, in_temp(:,iN), this.tOut);
    end
end
inputs(isnan(inputs)) = 0;

%             this.t_in = t_in;
%             inputs = in_temp;
this.inputs = inputs;
end




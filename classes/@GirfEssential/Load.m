function Load(this, filename)
%Function to load girf data from file into GirfEssential object
%
% IN
% filename  Name of file to load data from
%
% OUT
%
% EXAMPLE
%   girfE.Load(mySavedGirfFilename);
%
%   See also GirfEssential
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%


if ~exist(filename,'file')
    error('Could not find file: %s',filename)
else
    in = load(filename);
    loadedVars = fieldnames(in);
    classProps = properties(this);
    for iV = 1:length(classProps)
        if ismember(classProps{iV},loadedVars)
            mp = findprop(this,classProps{iV});
            if ~mp.Dependent
                this.(classProps{iV}) = in.(classProps{iV});
            end
        end
    end
end

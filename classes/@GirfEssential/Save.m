function Save(this, filename, overwrite)
%Function to save data in GirfEssential object to file
%
% IN
% filename  Name of file to save data into
% overwrite Overwrite existing file?
%
% OUT
%
% EXAMPLE
%   girfE.Save(mySaveGirfFilename);
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


if nargin < 3
    overwrite = 0;
end

if exist(filename,'file')
    if ~overwrite
        warning('Already existing file: %s, not saving girf!!',filename)
        return
    end
end

classProps = properties(this);
firstSave = 1;
for iV = 1:length(classProps)
    mp = findprop(this,classProps{iV});
    if ~mp.Dependent
        eval([classProps{iV} ' = this.' classProps{iV} ';']);
        if firstSave == 1
            save(filename,classProps{iV});
            firstSave = 0;
        else
            save(filename,classProps{iV},'-append');
        end
    end
end
end

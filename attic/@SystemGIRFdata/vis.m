function vis(this, varargin)
% Visualization function
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


% Assign input variables and defaults
argin_names = {'GIRFs' 'domain' 'k_plot' 'BW' 'BW_rc' 'integration'};
defaults = {'{''X'' ''Y'' ''Z''}' 'f' 'self' '0' '0' '0'};
for iN = 1:nargin
    eval([argin_names{iN} ' = varargin{' num2str(iN) '};']);
end
for iN = nargin+1:length(argin_names)
    eval([argin_names{iN} ' = ' defaults{iN} ';']);
end

switch domain % TODO: f and t in SystemGIRFdata
    case 'f'
        xdata = this.f;
        if BW
            inds = find(abs(xdata)<this.BW/2); % BW is full width
        else
            inds = 1:length(this.f);
        end
        xdata = xdata(inds);
    case 't'
        xdata = centered_time(this.dt,length(this.t_c));
        if integration
            xdata = xdata +this.dt/2;
        end
end

% Make figure and plot according to type
figure('Color', [1 1 1], 'Name', ['GIRFs, ' domain])
n_GIRFs = length(GIRFs);
for iG = 1:n_GIRFs
    if strcmp(k_plot,'self')
        k_GIRF = this.(GIRFs{iG}).self;
    end
    
    if strcmp(domain,'f')
        ax(1) = subplot(2,1,1);
        plot(xdata, abs(this.(GIRFs{iG}).GIRF(inds,k_GIRF)))
        hold all
        ax(2) = subplot(2,1,2);
        plot(xdata, centered_phase(this.(GIRFs{iG}).GIRF(inds,k_GIRF)))
        hold all
        linkaxes(ax, 'x')
    else
        if integration
            ydata = this.gamma*this.dt*cumsum(BW_filter(this.(GIRFs{iG}).girft(:,k_GIRF),xdata,BW,'rc',BW_rc));
        else
            ydata = BW_filter(this.(GIRFs{iG}).girft(:,k_GIRF),xdata,BW,'rc',BW_rc);
        end
        plot(xdata,ydata)
        hold all
    end
    
end
legend(GIRFs)
end
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
argin_names = {'plottype' 'domain' 'k_plot' 'BW' 'BW_rc' 'pulses' 'integration'};
defaults = {'''GIRF''' '''f''' '''self''' '0' '0' '''all''' '0'};
for iN = 1:nargin-1
    eval([argin_names{iN} ' = varargin{' num2str(iN) '};']);
end
for iN = nargin:length(argin_names)
    eval([argin_names{iN} ' = ' defaults{iN} ';']);
end
if strcmp(k_plot,'self')
    k_plot = this.self;
end
if strcmp(pulses, 'all')
    pulses = 1:this.nIn;
end

% Make vector for x-axis in plot
switch domain
    case 'f'
        xdata = this.f;
        if isempty(xdata)
            xdata = time2freq(this.tOut);
            this.f = xdata;
        end
        if BW
            inds = find(abs(xdata)<this.BW/2); % BW is full width
        else
            inds = 1:length(this.f);
        end
        xdata = xdata(inds);
    case 't'
        if integration
            xdata = this.t_k;
        else
            xdata = this.tOut;
        end
        
end

% Make figure and plot according to type
figure('Color', [1 1 1], 'Name', [plottype ' ' this.channel ', ' domain])
switch plottype
    case {'GIRF' 'girf'}
        if strcmp(domain,'f')
            ax(1) = subplot(2,1,1);
            plot(xdata, abs(this.GIRF(inds,k_plot)))
            ax(2) = subplot(2,1,2);
            plot(xdata, CenteredPhase(this.GIRF(inds,k_plot)))
            linkaxes(ax, 'x')
            if ~isempty(this.BW) && this.BW ~= 0
                set(ax,'XLim',[-this.BW/2 this.BW/2])
            end
        else
            xdata = CenteredTime(this.dtOut,size(this.GIRF,1));
            if integration
                xdata = xdata +this.dtOut/2;
                ydata = this.gamma*this.dtOut*cumsum(BWfilter(this.girft(:,k_plot),xdata,BW,'rc',BW_rc));
            else
                ydata = BWfilter(real(this.girft(:,k_plot)),xdata,BW,'rc',BW_rc);
            end
            plot(xdata,ydata)
        end
    case {'IN' 'in'}
        if strcmp(domain,'f')
            ax(1) = subplot(2,1,1);
            plot(xdata, abs(this.IN(inds,pulses)))
            ax(2) = subplot(2,1,2);
            plot(xdata, CenteredPhase(this.IN(inds,pulses)))
            linkaxes(ax, 'x')
        else
            if integration
                ydata = this.gamma*this.dt*cumsum(BWfilter(this.inputs(:,pulses),xdata,BW,'rc',BW_rc));
            else
                ydata = BWfilter(this.inputs(:,pulses),xdata,BW,'rc',BW_rc);
            end
            plot(xdata, ydata)
        end
    case {'OUT' 'out'}
        if strcmp(domain,'f')
            ax(1) = subplot(2,1,1);
            plot(xdata, abs(this.OUT(inds,pulses,k_plot)))
            ax(2) = subplot(2,1,2);
            plot(xdata, CenteredPhase(this.OUT(inds,pulses,k_plot)))
            linkaxes(ax, 'x')
        else
            if integration
                ydata = this.gamma*this.dt*cumsum(BWfilter(this.outputs(:,pulses,k_plot),xdata,BW,'rc',BW_rc));
            else
                ydata = BWfilter(this.outputs(:,pulses,k_plot),xdata,BW,'rc',BW_rc);
            end
            plot(xdata, ydata)
        end
    case {'INOUT' 'inout'}
        if strcmp(domain,'f')
            ax(1) = subplot(2,1,1);
            plot(xdata, abs(this.IN(inds,pulses)), 'k')
            hold all
            plot(xdata, abs(this.OUT(inds,pulses,k_plot)))
            ax(2) = subplot(2,1,2);
            plot(xdata, CenteredPhase(this.IN(inds,pulses)), 'k')
            hold all
            plot(xdata, CenteredPhase(this.OUT(inds,pulses,k_plot)))
            linkaxes(ax, 'x')
        else
            if integration
                ydata_in = this.gamma*this.dt*cumsum(BWfilter(this.inputs(:,pulses),xdata,BW,'rc',BW_rc));
                ydata_out = this.gamma*this.dt*cumsum(BWfilter(this.outputs(:,pulses,k_plot),xdata,BW,'rc',BW_rc));
            else
                ydata_in = BWfilter(this.inputs(:,pulses),xdata,BW,'rc',BW_rc);
                ydata_out = BWfilter(this.outputs(:,pulses,k_plot),xdata,BW,'rc',BW_rc);
            end
            plot(xdata, ydata_in, 'k')
            hold all
            plot(xdata, ydata_out)
        end
    case 'sensitivity'
        IN_w = zeros(size(this.inputs));
        if ~iscell(this.dyns) && ~isempty(this.dyns)
            IN_w = sqrt(length(this.dyns))*this.IN;
        elseif ~isempty(this.dyns)
            for iN = 1:length(this.nIn)
                IN_w(:,iN) = sqrt(length(this.dyns{iN}))*this.IN(:,iN);
            end
        else
            IN_w = this.IN;
        end
        IN_rss = sqrt(sum(abs(IN_w).^2,2));
        plot(xdata, abs(IN_w(inds,:)))
        hold all
        plot(xdata, abs(IN_rss(inds)), 'k')
end
end
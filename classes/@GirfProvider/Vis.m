function fH = Vis(this, plottype, domain, inChannel, outBasis, waveforms, integration)
% Visualizes various GIRF variables (GIRF, input, output, etc)
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
if nargin < 7 || isempty(integration)
    integration = 0;
end
if nargin < 6 || isempty(waveforms)
    waveforms = 1:size(this.in,3);
end
if nargin < 5 || isempty(outBasis)
    outBasis = this.selfBasis;
end
if nargin < 4 || isempty(inChannel)
    inChannels = this.inChannels;
end
if nargin < 3  || isempty(domain)
    domain = 'f';
end
if nargin < 2  || isempty(plottype)
    plottype = 'GIRF';
end

% Make figure and plot according to type
for iIn = 1:length(inChannels)
    fH = figure('Color', [1 1 1], 'Name', [plottype ' ' inChannels{iIn} ', ' domain]);
    switch plottype
        case {'GIRF' 'Girf' 'girf'} % Plot computed GIRF ...
            if strcmp(domain,'f') % ... in the frequency domain
                if isempty(this.girf) || isempty(this.freq)
                    this.ConvertDomain('t2f');
                end
                f = this.freq;
                girf = this.girf(:,outBasis(iIn),iIn);
                ax(1) = subplot(2,1,1);
                plot(f, abs(girf))
                ax(2) = subplot(2,1,2);
                plot(f, CenteredPhase(girf))
                linkaxes(ax, 'x')
            else  % ... in the time domain
                if isempty(this.girfTime) || isempty(this.time)
                    this.ConvertDomain('f2t');
                end
                time = this.time;
                girfTime = this.girfTime(:,outBasis(iIn),iIn);
                if integration
                    time = time +this.dt/2;
                    girfTime = this.dt*cumsum(girfTime);
                end
                plot(time,girfTime)
            end
        case {'IN' 'in'} % Plot input pulses ...
            if strcmp(domain,'f')  % ... in the frequency domain
                f = time2freq(this.timeIn);
                ax(1) = subplot(2,1,1);
                plot(f, abs(squeeze(this.inFreq(:,iIn,waveforms))))
                ax(2) = subplot(2,1,2);
                plot(f, CenteredPhase(squeeze(this.inFreq(:,iIn,waveforms))))
                linkaxes(ax, 'x')
            else  % ... in the time domain
                t = this.timeIn;
                in = squeeze(this.in(:,iIn,waveforms));
                if integration
                    t = t + this.dtIn/2;
                    in = this.dtIn*cumsum(in);
                end
                plot(t, in)
            end
        case {'OUT' 'out'} % Plot output pulses ...
            if strcmp(domain,'f') % ... in the frequency domain
                f = time2freq(this.timeOut);
                ax(1) = subplot(2,1,1);
                plot(f, abs(squeeze(this.outFreq(:,outBasis(iIn),waveforms))))
                ax(2) = subplot(2,1,2);
                plot(f, CenteredPhase(squeeze(this.outFreq(:,outBasis(iIn),waveforms))))
                linkaxes(ax, 'x')
            else % ... in the time domain
                t = this.timeOut;
                out = squeeze(this.out(:,outBasis(iIn),waveforms));
                if integration
                    t = t + this.dtOut/2;
                    out = this.dtOut*cumsum(out);
                end
                plot(t, out)
            end
        case {'INOUT' 'inout'} % Plot output pulses ...
            if strcmp(domain,'f') % ... in the frequency domain
                fIn = time2freq(this.timeIn);
                fOut = time2freq(this.timeOut);
                in = squeeze(this.inFreq(:,iIn,waveforms));
                out = squeeze(this.outFreq(:,outBasis(iIn),waveforms));
                ax(1) = subplot(2,1,1);
                plot(fIn, abs(in), 'k')
                hold all
                plot(fOut, abs(out))
                ax(2) = subplot(2,1,2);
                plot(fIn, CenteredPhase(in), 'k')
                hold all
                plot(fOut, CenteredPhase(out))
                linkaxes(ax, 'x')
            else % ... in the time domain
                tIn = this.timeIn;
                tOut = this.timeOut;
                in = squeeze(this.in(:,iIn,waveforms));
                out = squeeze(this.out(:,outBasis(iIn),waveforms));
                if integration
                    tIn = tIn + this.dtIn/2;
                    tOut = tOut + this.dtOut/2;
                    in = this.dtIn*cumsum(in);
                    out = this.dtOut*cumsum(out);
                end
                plot(tIn, in, 'k')
                hold all
                plot(tOut, out)
            end
        case 'sensitivity' % Plot sensitivity of GIRF measurement
            f = time2freq(this.timeIn);
            inFreq = squeeze(this.inFreq(:,iIn,waveforms));
            inFreqRSS = sqrt(sum(abs(inFreq).^2,2));
            plot(f, abs(inFreq))
            hold all
            plot(f, abs(inFreqRSS), 'k')
    end
end
end
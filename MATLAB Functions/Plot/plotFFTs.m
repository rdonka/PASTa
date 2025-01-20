function [allffts] = plotFFTs(data,whichfile,maintitle,whichfs,varargin)
% PLOTFFTs      Creates frequency magnitude plots of fiber photometry 
%               streams. This function will take the FFTs and plot the 
%               streams sig, baq, baq_scaled, sigsub, and sigfilt. Use this 
%               function in a loop to make plots for all sessions in the 
%               data structure.
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       streams to be plotted.
%
%       WHICHFILE:      The file number to plot (this can be set in a for
%                       loop to plot all files).
%
%       MAINTITLE:      The main title for the overall plot to display
%                       above the individual tiles. For example,
%                       '427 - Treatment: Morphine'
%
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
% OPTIONAL INPUTS:
%       XMAX:           Frequency cutoff for the x axis. Frequencies above
%                       this value will be excluded from the plots. To plot
%                       all frequencies, set to 'actual'. Default: 20.
%
%       SAVEOUTPUT:     Set to 1 to automatically save trace plots as png 
%                       to the plot file path. Default: 0.
%
%       PLOTFILEPATH:   Required if SAVEOUTPUT is set to 1. The specific 
%                       path to save the plot to. Note that this must be 
%                       the entire path from computer address to folder, 
%                       ending in the filename for the specific plot.
%                       For example: 
%                       'C:\Users\rmdon\Box\Injection Transients\Figures\SessionFFTs_427_Morphine.png'
%
% OUTPUT:
%       ALLTRACES:      A plot object containing subplots for each input
%                       stream.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

disp(append('PLOTFFTS: Plotting frequency magnitude plots for file: ',num2str(whichfile)))

%% Prep FFTs
[sigFFT,sigF] = preparestreamFFT(data(whichfile).sig,data(whichfile).(whichfs));
[baqFFT,baqF] = preparestreamFFT(data(whichfile).baq,data(whichfile).(whichfs));
[baqscaledFFT,baqscaledF] = preparestreamFFT(data(whichfile).baq_scaled,data(whichfile).(whichfs));
[sigsubFFT,sigsubF] = preparestreamFFT(data(whichfile).sigsub,data(whichfile).(whichfs));
[sigfiltFFT,sigfiltF] = preparestreamFFT(data(whichfile).sigfilt,data(whichfile).(whichfs));

%% Prepare Inputs
    inputs = struct(...
        'saveoutput',[],...
        'plotfilepath',[],...
        'xmax',[]);
    inputs = parseArgsLite(varargin,inputs);

    % Prepare defaults and check for optional inputs
    if isempty(inputs.saveoutput)
        saveoutput = 0; % Defaults to 0 - skip saving plots
        inputs.saveoutput = saveoutput;
    else
        saveoutput = inputs.saveoutput;
    end
    if isempty(inputs.xmax)
        xmax = 20; % Defaults to 20 - cuts off the x axis of the plot at 20 Hz
        inputs.xmax = xmax;
        disp(append('  X axis (frequency) max cut off at: ',num2str(xmax),' hz'))
    elseif inputs.xmax == 'actual'
        xmax = max(sigF);
        disp(append('  X axis (frequency) max set to actual max: ',num2str(xmax),' hz'))
    else
        xmax = inputs.xmax;
        disp(append('  X axis (frequency) max cut off at: ',num2str(xmax),' hz'))
    end

    if isempty(inputs.plotfilepath) & saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        disp(' ERROR: SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

%% Prep colors
    sigcolor = '#0092FF';
    baqcolor = '#8200C8';
    baqscaledcolor = '#e70089';
    sigsubcolor = '#00C296';
    sigfiltcolor = '#4CBB17';
    
%% Prep axis variables
    currxlength = sum(sigF<xmax);
    currxticksize = floor(xmax/20); % Find total number of minutes per session - helper variable to determine ticks
    currxticks = 0:currxticksize:xmax;

    fftymax = 10^1;
    fftymin = 10^-7;
    fftyticks = [10^-6, 10^-3, 10^0];

%% Plot traces
    close all
    ntraces = 6;

    % Create tiled layout
    allffts = tiledlayout(ntraces, 1, 'Padding','compact', 'TileSpacing','compact');

    
    % Plot raw signal
    nexttile;
    hold on;
    plot(sigF(1:currxlength), sigFFT(1:currxlength), 'Color', sigcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Raw Signal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;

    % Plot raw background
    nexttile;
    hold on;
    plot(baqF(1:currxlength), baqFFT(1:currxlength), 'Color', baqcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Raw Background');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;

    % Plot raw signal with raw background
    nexttile;
    hold on;
    plot(sigF(1:currxlength), sigFFT(1:currxlength), 'Color', sigcolor);
    plot(baqF(1:currxlength), baqFFT(1:currxlength), 'Color', baqcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Raw Signal with Raw Background');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;

    % Plot raw signal with scaled background
    nexttile;
    hold on;
    plot(sigF(1:currxlength), sigFFT(1:currxlength), 'Color', sigcolor);
    plot(baqscaledF(1:currxlength), baqscaledFFT(1:currxlength), 'Color', baqscaledcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Raw Signal with Scaled Background');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;

    % Plot subtracted signal
    nexttile;
    hold on;
    plot(sigsubF(1:currxlength), sigsubFFT(1:currxlength), 'Color', sigsubcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Subtracted Signal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;
    
    % Plot filtered signal
    nexttile;
    hold on;
    plot(sigfiltF(1:currxlength), sigfiltFFT(1:currxlength), 'Color', sigfiltcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Subtracted and Filtered Signal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;
    
    % Add a main title for the entire tiled layout
    title(allffts, maintitle);

    if saveoutput == 1
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntraces]);
        exportgraphics(gcf,append(plotfilepath, '.png'),'Resolution',300)
    end
end

% Copyright (C) 2024 Rachel Donka
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
function [alltraces] = plotTraces(data,whichfile,maintitle,varargin)
% PLOTTRACES    Plots whole session fiber photometry traces. This function
%               will plot the streams sig, baq, baq_scaled, sigsub, and
%               sigfilt. Use this function in a loop to plot streams for
%               all sessions in the data structure.
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
%                       '427 - Treatment: Morphine''
%
% OPTIONAL INPUTS:
%       SAVEOUTPUT:     Set to 1 to automatically save trace plots as png 
%                       to the plot file path. Default: 0.
%
%       PLOTFILEPATH:   Required if SAVEOUTPUT is set to 1. The specific 
%                       path to save the plot to. Note that this must be 
%                       the entire path from computer address to folder, 
%                       ending in the filename for the specific plot.
%                       For example: 
%                       'C:\Users\rmdon\Box\Injection Transients\Figures\SessionTraces_427_Morphine.png'
%
% OUTPUT:
%       ALLTRACES:      A plot object containing subplots for each input
%                       stream.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Inputs
    inputs = struct(...
        'saveoutput',[],...
        'plotfilepath',[]);
    inputs = parseArgsLite(varargin,inputs);

    disp(append('PLOTTRACES: Plotting session traces for file: ',num2str(whichfile)))

    % Prepare defaults and check for optional inputs
    if isempty(inputs.saveoutput)
        saveoutput = 0; % Defaults to 0 - skip saving plots
        inputs.saveoutput = saveoutput;
    else
        saveoutput = inputs.saveoutput;
    end

    if isempty(inputs.plotfilepath) & saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        disp('   ERROR: SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

%% Prep colors
    sigcolor = '#0092FF';
    baqcolor = '#8200C8';
    sigsubcolor = '#00C296';
    sigfiltcolor = '#4CBB17';
    
%% Prep axis variables
    currxlength = length(data(whichfile).sig);
    currxmins = (length(data(whichfile).sig)/data(whichfile).fs)/60; % Find total number of minutes per session - helper variable to determine ticks
    currxticklabels = 0:5:floor(currxmins/5)*5;
    currxticks = floor(currxticklabels.*60.*data(whichfile).fs); % Determine x axis ticks - add ticks every 5 minutes

    ymaxsig = ceil(max(data(whichfile).sig)+(0.1*max(data(whichfile).sig)));
    yminsig = floor(min(data(whichfile).sig)-(0.1*min(data(whichfile).sig)));
    yticksizesig = round((ymaxsig-yminsig)/4,0); % Find size of ticks to generate 5 y axis ticks total
    currytickssig = yminsig:yticksizesig:ymaxsig;

    ymaxbaq = ceil(max(data(whichfile).baq)+(0.1*max(data(whichfile).baq)));
    yminbaq = floor(min(data(whichfile).baq)-(0.1*min(data(whichfile).baq)));
    yticksizebaq = round((ymaxbaq-yminbaq)/4,0); % Find size of ticks to generate 5 y axis ticks total
    curryticksbaq = yminbaq:yticksizebaq:ymaxbaq;

    ymaxsigsub = ceil(max(data(whichfile).sigsub)+(0.1*max(data(whichfile).sigsub)));
    yminsigsub = floor(min(data(whichfile).sigsub)-(0.1*min(data(whichfile).sigsub)));
    yticksizesigsub = round((ymaxsigsub-yminsigsub)/4,1); % Find size of ticks to generate 5 y axis ticks total
    currytickssigsub = yminsigsub:yticksizesigsub:ymaxsigsub;

    ymaxsigfilt = ceil(max(data(whichfile).sigfilt)+(0.1*max(data(whichfile).sigfilt)));
    yminsigfilt = floor(min(data(whichfile).sigfilt)-(0.1*min(data(whichfile).sigfilt)));
    yticksizesigfilt = round((ymaxsigfilt-yminsigfilt)/4,1); % Find size of ticks to generate 5 y axis ticks total
    currytickssigfilt = yminsigfilt:yticksizesigfilt:ymaxsigfilt;

%% Plot traces
    close all
    ntraces = 5;

    % Create tiled layout
    alltraces = tiledlayout(ntraces, 1, 'Padding','compact', 'TileSpacing','compact');

    % Plot raw signal
    nexttile;
    hold on;
    plot(data(whichfile).sig, 'Color', sigcolor);
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([yminsig ymaxsig]);
    yticks(currytickssig);
    title('Raw Signal');
    xlabel('Minute');
    ylabel('Fluorescence (A.U.)');
    hold off;
    
    % Plot raw background
    nexttile;
    hold on;
    plot(data(whichfile).baq, 'Color', baqcolor);
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([yminbaq ymaxbaq]);
    yticks(curryticksbaq);
    title('Raw Background');
    xlabel('Minute');
    ylabel('Fluorescence (A.U.)');
    hold off;
    
    % Plot raw signal with scaled background
    nexttile;
    hold on;
    plot(data(whichfile).sig, 'Color', sigcolor);
    plot(data(whichfile).baq_scaled, 'Color', baqcolor);
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([yminsig ymaxsig]);
    yticks(currytickssig);
    title('Raw Signal with Scaled Background');
    xlabel('Minute');
    ylabel('Fluorescence (A.U.)');
    hold off;
    
    % Plot subtracted signal
    nexttile;
    hold on;
    plot(data(whichfile).sigsub, 'Color', sigsubcolor);
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([yminsigsub ymaxsigsub]);
    yticks(currytickssigsub);
    title('Subtracted Signal');
    xlabel('Minute');
    ylabel('dF/F');
    hold off;
    
    % Plot filtered signal
    nexttile;
    hold on;
    plot(data(whichfile).sigfilt, 'Color', sigfiltcolor);
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([yminsigfilt ymaxsigfilt]);
    yticks(currytickssigfilt);
    title('Subtracted and Filtered Signal');
    xlabel('Minute');
    ylabel('dF/F');
    hold off;
    
    % Add a main title for the entire tiled layout
    title(alltraces, maintitle);
    fontsize(alltraces,scale=.8)

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
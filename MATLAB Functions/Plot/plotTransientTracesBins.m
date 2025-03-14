function [alltransienttracebins] = plotTransientTracesBins(data,whichfile,maintitle,whichfs,whichtransients,varargin)
% PLOTTRANSIENTTRACESBINS    Plots traces of each transient by bin.
%
% Copyright (C) 2025 Rachel Donka. Licensed under the GNU General Public License v3.
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
%       WHICHTRANSIENTS:    String; The name of the parent field containing 
%                           the table of transients that you want to identify 
%                           bins for. For example, 'sessiontransients_blmin_3SD'.
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
%       ALLTRANSIENTTRACEBINS: A plot object containing subplots for each
%                              bin.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Inputs
    inputs = struct(...
        'saveoutput',[],...
        'plotfilepath',[]);
    inputs = parseArgsLite(varargin,inputs);

    disp(append('PLOTTRANSIENTTRACESBINS: Plotting session traces for file: ',num2str(whichfile)))

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
    nbins = data(whichfile).(whichtransients).BinSettings.nbins;
    binfieldname = append('Bin_',num2str(data(whichfile).(whichtransients).BinSettings.binlengthmins));
    
    currxlength = length(data(whichfile).(whichtransients).transientstreamdata);
    % currxseconds = currxlength/data(whichfile).(whichfs); % Find total number of minutes per session - helper variable to determine ticks
    currxticklabels = 0:.5:currxseconds;
    currxticks = currxticklabels.*data(whichfile).(whichfs); % Determine x axis ticks - add ticks every 5 minutes

    ymax = ceil(max(data(whichfile).(whichtransients).transientstreamdata, [], 'all')+(0.1*max(data(whichfile).(whichtransients).transientstreamdata, [], 'all')));
    ymin = floor(min(data(whichfile).(whichtransients).transientstreamdata, [], 'all')-(0.1*min(data(whichfile).(whichtransients).transientstreamdata, [], 'all')));
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    curryticks = ymin:yticksize:ymax;

%% Plot traces
    close all

    % Create tiled layout
    alltransienttracebins = tiledlayout(ceil(nbins/5), 5, 'Padding','compact', 'TileSpacing','compact');

    % Plot raw signal
    for eachbin = 1:nbins
        bintransients = find(data(whichfile).(whichtransients).transientquantification.(binfieldname) == eachbin);
        nexttile;
        set(gca, 'ColorOrder', lines(length(bintransients)), 'NextPlot', 'replacechildren');
        hold on;
        plot(data(whichfile).(whichtransients).transientstreamdata(bintransients,:)');
        xlim([0 currxlength]);
        xticks(currxticks);
        xticklabels(currxticklabels);
        ylim([ymin ymax]);
        yticks(curryticks);
        title(append('Bin ', num2str(eachbin)));
        xlabel('Seconds');
        ylabel('Z Score');
        hold off;
    end
   
    % Add a main title for the entire tiled layout
    title(alltransienttracebins, maintitle, 'Interpreter', 'none');

    if saveoutput == 1
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 6, 1.75*ntraces]);
        exportgraphics(gcf,append(plotfilepath, '.png'),'Resolution',300)
    end
end

% Copyright (C) 2025 Rachel Donka
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
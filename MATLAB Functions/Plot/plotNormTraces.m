function [normtraces] = plotNormTraces(data,whichfile,whichstreams,whichfs,maintitle,streamtitles,varargin)
% PLOTTRACES    Plots whole session fiber photometry traces. This function
%               will plot the streams sig, baq, baq_scaled, sigsub, and
%               sigfilt. Use this function in a loop to plot streams for
%               all sessions in the data structure.
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
%       WHICHSTREAMS:   A cell array containing the names (string) of the 
%                       fields in the data structure containing the normalized 
%                       (Z scored) stream to be plotted.
%
%       WHICHFS:        String; The name of the field that contains the
%                       sampling rate of the stream used for transient
%                       detection. For example, 'fs'.
%
%       MAINTITLE:      The main title for the overall plot to display
%                       above the individual tiles. For example,
%                       '427 - Treatment: Morphine'
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
%       NORMTRACE:      A plot object containing subplots for each input
%                       stream.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Inputs
    inputs = struct(...
        'saveoutput',[],...
        'plotfilepath',[]);
    inputs = parseArgsLite(varargin,inputs);

    disp(append('PLOTNORMTRACES: Plotting normalized traces for file: ',num2str(whichfile)))
    disp(append('   whichstreams: '))
    disp(whichstreams)

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
    signormcolor = '#1300b6';
    
    %% Prep tiled layout
    ntiles = length(whichstreams); % Number of tiles
    
    close all
    normtraces = tiledlayout(ntiles, 1, 'Padding','compact', 'TileSpacing','compact');
        
    %% Prep y axis variables
    streammax = 0;
    streammin = 0;
    for eachstream = 1:length(whichstreams)
        currstream = char(whichstreams(eachstream));
        % Find max
        currstreammax = max(data(whichfile).(currstream));
        if currstreammax > streammax
            streammax = currstreammax;
        end
        % Find min
        currstreammin = min(data(whichfile).(currstream));
        if currstreammin < streammin
            streammin = currstreammin;
        end
    end
    
    ymax = ceil(max(data(whichfile).(currstream))+(0.1*max(data(whichfile).(currstream))));
    ymin = floor(streammin-(0.1*streammin));
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    yticklocs = ymin:yticksize:ymax;

    %% Plot each stream
    for eachstream = 1:length(whichstreams)
        currstream = char(whichstreams(eachstream));
        currstreamtitle = char(streamtitles(eachstream));
        
        % Prep x axis variables
        currxlength = length(data(whichfile).(currstream));
        currxmins = (length(data(whichfile).(currstream))/data(whichfile).(whichfs))/60; % Find total number of minutes per session - helper variable to determine ticks
        currxticklabels = 0:5:floor(currxmins/5)*5;
        currxticks = floor(currxticklabels.*60.*data(whichfile).(whichfs)); % Determine x axis ticks - add ticks every 5 minutes
    
        % Plot normalized signal
        nexttile;
        hold on;
        plot(data(whichfile).(currstream), 'Color', signormcolor);
        xlim([0 currxlength]);
        xticks(currxticks);
        xticklabels(currxticklabels);
        ylim([ymin ymax]);
        yticks(yticklocs);
        title(currstreamtitle);
        xlabel('Minute');
        ylabel('Z Score');
        hold off;
    end

    % Add a main title for the entire tiled layout
    title(normtraces, maintitle);

    if saveoutput == 1
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntiles]);
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
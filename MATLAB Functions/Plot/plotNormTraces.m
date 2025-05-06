function [normtraces] = plotNormTraces(data,fileindex,streamfieldnames,fsfieldname,maintitle,streamtitles,varargin)
% PLOTNORMTRACES  Plot normalized (Z-scored) fiber photometry traces for 
%                 a session.
%
%   PLOTNORMTRACES(DATA, FILEINDEX, STREAMFIELDNAMES, FSFIELDNAME, MAINTITLE, STREAMTITLES, 'PARAM1', VAL1, ...)
%   plots the specified normalized data streams for a given session. This 
%   function can be used within a loop to plot streams for all sessions in 
%   the data structure.
%
% REQUIRED INPUTS:
%   DATA              - Structure array; must contain the streams to be plotted.
%
%   FILEINDEX         - Integer; index of the file (session) to plot. 
%                       This can be set within a loop to plot all files.
%
%   STREAMFIELDNAMES  - Cell array of strings; names of the fields in the 
%                       data structure containing the normalized (Z-scored) 
%                       streams to be plotted.
% 
%   FSFIELDNAME       - String; name of the field that contains the sampling 
%                       rate of the stream used for transient detection. 
%                       For example, 'fs'.
%
%   MAINTITLE         - String; main title for the overall plot, displayed 
%                       above the individual subplots. For example, 
%                       '427 - Treatment: Morphine'.
%
%   STREAMTITLES      - Cell array of strings; titles for each subplot 
%                       corresponding to each stream in STREAMFIELDNAMES.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'saveoutput'    - Logical; set to true to automatically save trace 
%                     plots to the specified plot file path.
%                     Default: false.
%

%   'plotfilepath'  - String; required if 'saveoutput' is set to true. 
%                     Specifies the full path, including the filename,
%                     where the plot should be saved. For example:
%                     'C:\Users\rmdon\Box\Injection Transients\Figures\SessionTraces_427_Morphine'.
%
% OUTPUTS:
%   NORMTRACES      - Figure handle; handle to the figure containing 
%                     subplots for each input stream.
%
% EXAMPLE USAGE:
%   % Define the streams to plot and their titles:
%   streams = {'sigfiltz_normsession', 'sigfiltz_normbaseline'};
%   titles = {'Whole Session Normalized', 'Baseline Normalized'};
%
%   % Plot normalized traces for the first session without saving:
%   plotNormTraces(data, 1, streams, 'fs', 'Subject 427 - Baseline Session', titles);
%
%   % Plot and save normalized traces for the second session:
%   plotNormTraces(data, 2, streams, 'fs', 'Subject 427 - Treatment Session', titles, ...
%       'saveoutput', true, 'plotfilepath', 'C:\Plots\Session2');
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

 %% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'saveoutput', defaultparameters.saveoutput, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1]))); % saveoutput: input must be logical or numeric (either 0 or 1); set to 1 to save plot automatically
    addParameter(p, 'plotfilepath', defaultparameters.plotfilepath, @(x) ischar(x) || isstring(x)); % plotfilepath: defaults to empty unless input is specified

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Display
    disp(['PLOTNORMTRACES: Plotting normalized traces for file: ',num2str(fileindex)])
    disp('   streamfieldnames: ')
    disp(streamfieldnames)

    % Prepare defaults and check for optional inputs
    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

    %% Prep colors
    signormcolor = '#1300b6';
    
    %% Prep tiled layout
    ntiles = length(streamfieldnames); % Number of tiles
    
    close all
    normtraces = tiledlayout(ntiles, 1, 'Padding','compact', 'TileSpacing','compact');
        
    %% Prep y axis variables
    streammax = 0;
    streammin = 0;
    for eachstream = 1:length(streamfieldnames)
        currstream = char(streamfieldnames(eachstream));
        % Find max
        currstreammax = max(data(fileindex).(currstream));
        if currstreammax > streammax
            streammax = currstreammax;
        end
        % Find min
        currstreammin = min(data(fileindex).(currstream));
        if currstreammin < streammin
            streammin = currstreammin;
        end
    end
    
    ymax = ceil(max(data(fileindex).(currstream))+(0.1*max(data(fileindex).(currstream))));
    ymin = floor(streammin-(0.1*streammin));
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    yticklocs = ymin:yticksize:ymax;

    %% Plot each stream
    for eachstream = 1:length(streamfieldnames)
        currstream = char(streamfieldnames(eachstream));
        currstreamtitle = char(streamtitles(eachstream));
        
        % Prep x axis variables
        currxlength = length(data(fileindex).(currstream));
        currxmins = (length(data(fileindex).(currstream))/data(fileindex).(fsfieldname))/60; % Find total number of minutes per session - helper variable to determine ticks
        if currxmins > 10.01
            currxticklabels = 0:5:floor(currxmins/5)*5;
            currxticks = floor(currxticklabels.*60.*data(fileindex).fs); % Determine x axis ticks - add ticks every 5 minutes
        else
            currxticklabels = 0:1:currxmins;
            currxticks = floor(currxticklabels.*60.*data(fileindex).fs); % Determine x axis ticks - add ticks every 5 minutes
        end
    
        % Plot normalized signal
        nexttile;
        hold on;
        plot(data(fileindex).(currstream), 'Color', signormcolor);
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
    title(normtraces, maintitle, 'Interpreter', 'none');

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntiles]);
        exportgraphics(gcf,append(params.plotfilepath, '.',params.outputfiletype),'Resolution',300)
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
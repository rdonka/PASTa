function [trialtraces] = plotTrialTraces(data,fileindex,trialstreamfieldname,maintitle,varargin)

% PLOTTRANSIENTBINS     Plots every trial for a session with detected
%                       transient peaks marked by circles.
%
% REQUIRED INPUTS:
%   DATA            - Structure array; must contain the stream to be plotted.
%
%   FILEINDEX       - Integer; index of the file (session) to plot. This can be
%                     set within a loop to plot all files.
%
%   STREAMFIELDNAME - String; The name of the field containing the stream
%                     used for transient detection. For example, 'sigfiltz_normsession'.
%
%   TRANSIENTDATA   - Structure array of the output from FINDTRANSIENTS
%                     with the field 'transientquantification'.
%
%   BINFIELDNAME    - String; The name of the field in TRANSIENTDATA under
%                     the 'transientquantification' table that contains the
%                     trial IDs for each transient event. For example, 
%                     'Bin_5min'.
%
%   MAINTITLE       - String; main title for the overall plot, displayed above
%                     the individual subplots. For example, '427 - Treatment:
%                     Morphine'.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'saveoutput'    - Logical; set to true to automatically save trace plots
%                     to the specified plot file path. Default: false.
%
%   'outputfiletype'- String; File type extension to save the figure as.
%                     Options supported: 'png', 'jpg', 'tiff', 'eps', and 'pdf'.
%                     Default: 'png'.
%
%   'plotfilepath'  - String; required if 'saveoutput' is set to true.
%                     Specifies the full path, including the filename, where
%                     the plot should be saved. For example:
%                     'C:\Users\rmdon\Box\Injection Transients\Figures\TransientBinTraces_427_Morphine'.
%
% OUTPUT:
%       ALLBINS:      A plot object containing subplots for each input
%                       stream.
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
    addParameter(p, 'outputfiletype', defaultparameters.outputfiletype, @(x) ischar(x) && ismember(x, {'png', 'jpg', 'tiff', 'eps', 'pdf'})); % outputfiletype: file type to save plot as if saveoutput is set to 1
    addParameter(p, 'plotfilepath', defaultparameters.plotfilepath, @(x) ischar(x) || isstring(x)); % plotfilepath: defaults to empty unless input is specified

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;
    
    % Display
    disp(['PLOTTRIALTRACES: Plotting trial traces of ',trialstreamfieldname,' for file: ',num2str(fileindex)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

    %% Prep colors
    tracecolor = '#4CBB17';
    
    %% Prep axis variables
    fs = data(fileindex).fs;
    ntrials = height(data(fileindex).(trialstreamfieldname));
    
    Xmax = max(data(fileindex).(trialstreamfieldname)(:));     % max over all elements
    Xmin = min(data(fileindex).(trialstreamfieldname)(:));     % min over all elements

    trialymax = ceil(Xmax+(0.1*Xmax));
    trialymin = floor(Xmin-(0.1*Xmin));    
    trialyticksize = round((trialymax-trialymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    trialyticks = trialymin:trialyticksize:trialymax;

    if contains(trialstreamfieldname, 'z') == 1
        currylabel = 'Z Score';
    else
        currylabel = 'dF/F';
    end

    %% Plot trial traces with transients
    close all

    % Create tiled layout
    trialtraces = tiledlayout(ceil(ntrials/5), 5, 'Padding','compact', 'TileSpacing','compact');

     for trial = 0:(ntrials-1)
        trialidx = trial+1;

        nexttile
        hold on
        
        % Prepare x axis variables for stream trace
        currxlength = length(data(fileindex).(trialstreamfieldname)(trialidx,:));
        currxseconds = currxlength/fs; % Find total number of minutes per session - helper variable to determine ticks
        currxticklabels = 0:2:currxseconds;
        currxticks = currxticklabels.*fs; % Determine x axis ticks - add ticks every 5 minutes
    
        % Plot stream trace with detected transients
        plot(data(fileindex).(trialstreamfieldname)(trialidx,:), 'Color', tracecolor)
        axis([0 currxlength trialymin trialymax])
        xticks(currxticks)
        xticklabels(currxticklabels)
        yticks(trialyticks)
        xlabel('Second')
        ylabel(currylabel)
        title(append('Trial ',num2str(trial+1)))
        hold off
     end
    
    % Add a main title for the entire tiled layout
    title(trialtraces, maintitle, 'Interpreter', 'none');

    % Set axis variables
    allAxes = findall(gcf,'Type','axes');
    set(allAxes, 'TickDir', 'out')
    set(allAxes, 'XColor', '#000000', 'YColor', '#000000'); % Make axis lines black
    set(gcf,'Color', 'w')

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 10, 2*ceil(ntrials/5)]);
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
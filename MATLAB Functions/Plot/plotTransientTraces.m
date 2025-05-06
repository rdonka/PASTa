function [alltransienttraces] = plotTransientTraces(transientdata,fileindex,maintitle,varargin)
% PLOTTRANSIENTTRACES   Plots overlaid traces for every transient event in
%                       a session.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA   - Structure array of the output from FINDTRANSIENTS
%                     with the field 'transientstreamdata'.
%
%   FILEINDEX       - Integer; index of the file (session) to plot. This can be
%                     set within a loop to plot all files.
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
    disp(['PLOTTRANSIENTTRACES: Plotting overlaid transient traces for file: ',num2str(fileindex)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end


    %% Prep axis variables 
    fs = transientdata(fileindex).params.findTransients.fs;

    currxlength = length(transientdata(fileindex).transientstreamdata);
    currxseconds = ceil(currxlength/fs); % Find total number of minutes per session - helper variable to determine ticks
    currxticklabels = 0:1:currxseconds;
    currxticks = currxticklabels.*fs; % Determine x axis ticks - add ticks every 5 minutes

    ymax = ceil(max(transientdata(fileindex).transientstreamdata, [], 'all')+(0.1*max(transientdata(fileindex).transientstreamdata, [], 'all')));
    ymin = floor(min(transientdata(fileindex).transientstreamdata, [], 'all')-(0.1*min(transientdata(fileindex).transientstreamdata, [], 'all')));
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    curryticks = ymin:yticksize:ymax;
    
    if contains(transientdata(fileindex).params.findTransients.streamfield, 'z') == 1
        currylabel = 'Z Score';
    else
        currylabel = 'dF/F';
    end

    %% Plot traces
    close all

    % Create tiled layout
    alltransienttraces = tiledlayout(1, 1, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    set(gca, 'ColorOrder', lines(height(transientdata(fileindex).transientstreamdata)), 'NextPlot', 'replacechildren');
    hold on;
    plot(transientdata(fileindex).transientstreamdata');
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([ymin ymax]);
    yticks(curryticks);
    xlabel('Seconds');
    ylabel(currylabel);
    hold off;
   
    % Add a main title for the entire tiled layout
    title(alltransienttraces, maintitle, 'Interpreter', 'none');

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 6]);
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
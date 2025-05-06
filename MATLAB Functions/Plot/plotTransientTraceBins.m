function [alltransienttracebins] = plotTransientTraceBins(transientdata,fileindex,binfieldname,maintitle,varargin)

% PLOTTRANSIENTTRACEBINS   Plots overlaid traces for every transient event
%                          by bin for a session.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA   - Structure array of the output from FINDTRANSIENTS
%                     with the field 'transientquantification'.
%
%   FILEINDEX       - Integer; index of the file (session) to plot. This can be
%                     set within a loop to plot all files.
%
%   BINFIELDNAME    - String; The name of the field in TRANSIENTDATA under
%                     the 'transientquantification' table that contains the
%                     bin IDs for each transient event. For example, 
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
    disp(['PLOTTRANSIENTTRACEBINS: Plotting overlaid transient traces by bin for file: ',num2str(fileindex)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end


    %% Prep colors
    sigcolor = '#0092FF';
    baqcolor = '#8200C8';
    sigsubcolor = '#00C296';
    sigfiltcolor = '#4CBB17';
    
    %% Prep axis variables 
    fs = transientdata(fileindex).params.findTransients.fs;

    binsamples = transientdata(fileindex).params.binTransients.(binfieldname).binlengthsamples;
    binlengthmins = transientdata(fileindex).params.binTransients.(binfieldname).binlengthmins;
    nbins = transientdata(fileindex).params.binTransients.(binfieldname).nbins;
    
    ymax = ceil(max(transientdata(fileindex).transientstreamdata, [], 'all')+(0.1*max(transientdata(fileindex).transientstreamdata, [], 'all')));
    ymin = floor(min(transientdata(fileindex).transientstreamdata, [], 'all')-(0.1*min(transientdata(fileindex).transientstreamdata, [], 'all')));    
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    curryticks = ymin:yticksize:ymax;

    currxlength = length(transientdata(fileindex).transientstreamdata);
    currxseconds = currxlength/fs; % Find total number of minutes per session - helper variable to determine ticks
    currxticklabels = 0:2:currxseconds;
    currxticks = currxticklabels.*fs; % Determine x axis ticks - add ticks every 5 minutes
    
    if contains(transientdata(fileindex).params.findTransients.streamfield, 'z') == 1
        currylabel = 'Z Score';
    else
        currylabel = 'dF/F';
    end

    %% Plot traces
    close all

    % Create tiled layout
    alltransienttracebins = tiledlayout(ceil(nbins/5), 5, 'Padding','compact', 'TileSpacing','compact');

    % Plot raw signal
    for eachbin = 1:nbins
        bintransients = find(transientdata(fileindex).transientquantification.(binfieldname) == eachbin);
        nexttile;
        set(gca, 'ColorOrder', lines(length(bintransients)), 'NextPlot', 'replacechildren');
        hold on;
        plot(transientdata(fileindex).transientstreamdata(bintransients,:)');
        xlim([0 currxlength]);
        xticks(currxticks);
        xticklabels(currxticklabels);
        ylim([ymin ymax]);
        yticks(curryticks);
        title(append('Bin ', num2str(eachbin)));
        xlabel('Seconds');
        ylabel(currylabel);
        hold off;
    end

    % Add a main title for the entire tiled layout
    title(alltransienttracebins, maintitle, 'Interpreter', 'none');

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 10, 2*ceil(nbins/5)]);
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
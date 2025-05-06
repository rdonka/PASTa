function [alltransients] = plotTransients(data,fileindex,streamfieldname,fsfieldname,transientdata,maintitle,varargin)
% PLOTTRANSIENTS        Plots whole session data stream trace with detected
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
%   FSFIELDNAME       - String; name of the field that contains the sampling 
%                       rate of the stream used for transient detection. 
%                       For example, 'fs'.
%
%   TRANSIENTDATA   - Structure array of the output from FINDTRANSIENTS
%                     with the field 'transientquantification'.
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
    disp(['PLOTTRANSIENTS: Plotting whole session trace with detected transients for file: ',num2str(fileindex)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

    %% Prep colors
    tracecolor = '#4CBB17';
    transientcolor = '#0004ff';
    
    %% Prep variables
    pklocs = transientdata(fileindex).transientquantification.maxloc;

    currxlength = length(data(fileindex).(streamfieldname));
    currxmins = (length(data(fileindex).(streamfieldname))/data(fileindex).(fsfieldname))/60; % Find total number of minutes per session - helper variable to determine ticks
    if currxmins > 10.01
        currxticklabels = 0:5:floor(currxmins/5)*5;
        currxticks = floor(currxticklabels.*60.*data(fileindex).fs); % Determine x axis ticks - add ticks every 5 minutes
    else
        currxticklabels = 0:1:currxmins;
        currxticks = floor(currxticklabels.*60.*data(fileindex).fs); % Determine x axis ticks - add ticks every 5 minutes
    end

    ymax = ceil(max(data(fileindex).(streamfieldname))+(0.1*max(data(fileindex).(streamfieldname))));
    ymin = floor(min(data(fileindex).(streamfieldname))+(0.1*min(data(fileindex).(streamfieldname))));
    yticksize = round((ymax-ymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    yticklocs = ymin:yticksize:ymax;

    if contains(streamfieldname, 'z') == 1
        currylabel = 'Z Score';
    else
        currylabel = '\Delta F/F';
    end

    %% Plot trace with transients
    close all

    figure()
    hold on;
    plot(data(fileindex).(streamfieldname), 'Color', tracecolor);
    plot(pklocs,data(fileindex).(streamfieldname)(pklocs),'o', 'Color', transientcolor)
    xlim([0 currxlength]);
    xticks(currxticks);
    xticklabels(currxticklabels);
    ylim([ymin ymax]);
    yticks(yticklocs);
    xlabel('Minute');
    ylabel(currylabel);
    title(maintitle, 'Interpreter', 'none');
    hold off;

    alltransients = gcf;

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 15, 3]);
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
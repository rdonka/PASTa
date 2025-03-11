function [alltraces] = plotTraces(data,whichfile,maintitle,varargin)
% PLOTTRACES  Plot whole-session fiber photometry traces.
%
%   PLOTTRACES(DATA, WHICHFILE, MAINTITLE, 'PARAM1', VAL1, ...) plots the
%   specified data streams for a given session, including 'sig', 'baq',
%   'baqscaled', 'sigsub', and 'sigfilt'. This function can be used within
%   a loop to plot streams for all sessions in the data structure.
%
% REQUIRED INPUTS:
%   DATA        - Structure array; must contain the streams to be plotted.
%
%   WHICHFILE   - Integer; index of the file (session) to plot. This can be
%                 set within a loop to plot all files.
%
%   MAINTITLE   - String; main title for the overall plot, displayed above
%                 the individual subplots. For example, '427 - Treatment:
%                 Morphine'.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'saveoutput'    - Logical; set to true to automatically save trace plots
%                     as PNG files to the specified plot file path.
%                     Default: false.
%
%   'plotfilepath'  - String; required if 'saveoutput' is set to true.
%                     Specifies the full path, including the filename, where
%                     the plot should be saved. For example:
%                     'C:\Users\rmdon\Box\Injection Transients\Figures\SessionTraces_427_Morphine.png'.
%
% OUTPUTS:
%   ALLTRACES   - Figure handle; handle to the figure containing subplots
%                 for each input stream.
%
% EXAMPLE USAGE:
%   % Plot traces for the first session without saving:
%   plotTraces(data, 1, 'Subject 427 - Baseline Session');
%
%   % Plot and save traces for the second session:
%   plotTraces(data, 2, 'Subject 427 - Treatment Session', 'saveoutput', true, 'plotfilepath', 'C:\Plots\Session2.png');
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
    disp(['PLOTTRACES: Plotting session traces for file: ',num2str(whichfile)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
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
    plot(data(whichfile).baqscaled, 'Color', baqcolor);
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
    title(alltraces, maintitle, 'Interpreter', 'none');

    if params.saveoutput == 1
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntraces]);
        exportgraphics(gcf,append(params.plotfilepath, '.png'),'Resolution',300)
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
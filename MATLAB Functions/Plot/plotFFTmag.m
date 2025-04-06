function [allffts] = plotFFTmag(data,fileindex,maintitle,fsfieldname,varargin)
% PLOTFFTMAG    Creates frequency magnitude plots of fiber photometry streams.
%               This function computes and plots the Fast Fourier Transform (FFT)
%               of specified data streams, including 'sig', 'baq', 'baqscaled',
%               'sigsub', and 'sigfilt'. It is designed to be used in a loop to
%               generate plots for all sessions in the data structure.
%
%   [ALLFFTS] = PLOTFFTMAG(DATA, FILEINDEX, MAINTITLE, FSFIELDNAME, 'PARAM1', VAL1, ...)
%   computes and plots the FFT magnitude spectra for the specified data streams in
%   the given session.
%
% REQUIRED INPUTS:
%   DATA        - Structure array; contains the data streams to be plotted.
%
%   FILEINDEX   - Numeric; index specifying the session (file) to plot.
%
%   MAINTITLE   - String; main title for the overall plot, displayed above 
%                 the individual subplots. (e.g., '427 - Treatment: Morphine').
%
%   FSFIELDNAME - String; name of the field containing the sampling rate of
%                 the streams (e.g., 'fs').
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%   'xmax'          - Numeric; frequency cutoff for the x-axis. Frequencies 
%                     above this value will be excluded from the plots. To 
%                     plot all frequencies, set to 'actual'. Default: 100.
%
%   'saveoutput'    - Logical; set to true to automatically save trace plots 
%                     to the specified plot file path. Default: false.
%
%   'outputfiletype'- String; File type extension to save the figure as.
%                     Options supported: 'png', 'jpg', 'tiff', 'eps', and 'pdf'.
%                     Default: 'png'.
%
%   'plotfilepath'  - String; required if 'saveoutput' is set to true. 
%                     Specifies the full path, including the filename, 
%                     where the plot should be saved (e.g.,
%                       'C:\Users\rmdon\Box\Injection Transients\Figures\SessionFFTmag_427_Morphine').
%
% OUTPUT
%   ALLFFTS     - Figure object; contains subplots for each input stream's
%                 FFT power spectrum.
%
% EXAMPLE USAGE:
%   % Plot FFT magnitude spectra for session 1 with a frequency cutoff of 50 Hz
%   plotFFTmag(data, 1, 'Subject 427 - Session 1', 'fs', 'xmax', 50);
%
%   See also: plotFFTpower
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
    addParameter(p, 'xmax', 100,@(x) (isnumeric(x) && isscalar(x) && (x > 0)) || (ischar(x) && strcmp(x, 'actual'))); % xmax: Must either be numeric and greater than 0, or set to 'actual'
    addParameter(p, 'saveoutput', defaultparameters.saveoutput, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1]))); % saveoutput: input must be logical or numeric (either 0 or 1); set to 1 to save plot automatically
    addParameter(p, 'outputfiletype', defaultparameters.outputfiletype, @(x) ischar(x) && ismember(x, {'png', 'jpg', 'tiff', 'eps', 'pdf'})); % outputfiletype: file type to save plot as if saveoutput is set to 1
    addParameter(p, 'plotfilepath', defaultparameters.plotfilepath, @(x) ischar(x) || isstring(x)); % plotfilepath: defaults to empty unless input is specified

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;
    
    % Display
    disp(['PLOTFFTMAG: Plotting FFT magnitude plots for file: ',num2str(fileindex)])

    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

    %% Prep FFTs
    [sigFFT,sigF] = preparestreamFFT(data(fileindex).sig,data(fileindex).(fsfieldname));
    [baqFFT,baqF] = preparestreamFFT(data(fileindex).baq,data(fileindex).(fsfieldname));
    [baqscaledFFT,baqscaledF] = preparestreamFFT(data(fileindex).baqscaled,data(fileindex).(fsfieldname));
    [sigsubFFT,sigsubF] = preparestreamFFT(data(fileindex).sigsub,data(fileindex).(fsfieldname));
    [sigfiltFFT,sigfiltF] = preparestreamFFT(data(fileindex).sigfilt,data(fileindex).(fsfieldname));

%% Prep colors
    sigcolor = '#0092FF';
    baqcolor = '#8200C8';
    baqscaledcolor = '#e70089';
    sigsubcolor = '#00C296';
    sigfiltcolor = '#4CBB17';
    
%% Prep axis variables
    currxlength = sum(sigF<params.xmax);
    currxticksize = floor(params.xmax/20); % Find total number of minutes per session - helper variable to determine ticks
    currxticks = 0:currxticksize:params.xmax;

    fftymax = 10^2;
    fftymin = 10^-8;
    fftyticks = [10^-8, 10^-4, 10^0];

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
    xlim([-0.1 params.xmax]);
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
    xlim([-0.1 params.xmax]);
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
    xlim([-0.1 params.xmax]);
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
    xlim([-0.1 params.xmax]);
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
    xlim([-0.1 params.xmax]);
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
    xlim([-0.1 params.xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title('Subtracted and Filtered Signal');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (Log)');
    set(gca,'TickDir','out')
    hold off;
    
    % Add a main title for the entire tiled layout
    title(allffts, maintitle, 'Interpreter', 'none');

    if params.saveoutput == 1
        disp(['   Automatically saved as ', params.outputfiletype, ' to: ', params.plotfilepath,'.',params.outputfiletype])
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntraces]);
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
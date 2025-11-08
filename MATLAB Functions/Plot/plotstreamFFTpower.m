function [FFTplot] = plotstreamFFTpower(data,fileindex,maintitle,fsfieldname,streamfieldname,varargin)
% PLOTSTREAMFFTPOWER  Creates frequency power plots of the specified fiber
% photometry stream for one session.
%               This function computes and plots the Fast Fourier Transform (FFT)
%               of the specified data stream in STREAMFIELDNAME. It is designed to be used in a loop to
%               generate plots for all sessions in the data structure.
%
%   [FFTplot] = PLOTSTREAMFFTPOWER(DATA, FILEINDEX, MAINTITLE, FSFIELDNAME, STREAMFIELDNAME)
%   computes and plots the FFT power spectra for the specified data stream in
%   the given session.
%
% REQUIRED INPUTS:
%   DATA        - Structure array; contains the data stream to be plotted.
%
%   FILEINDEX   - Numeric; index specifying the session (file) to plot.
%
%   MAINTITLE   - String; main title for the overall plot, displayed above 
%                 the individual subplots. (e.g., '427 - Treatment: Morphine').
%
%   FSFIELDNAME - String; name of the field containing the sampling rate of
%                 the streams (e.g., 'fs').
%
%   STREAMFIELDNAME - String; name of field containing the stream to plot
%                   FFT power spectrum.
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
%                       'C:\Users\rmdon\Box\Injection Transients\Figures\SessionFFTpower_427_Morphine').
%
% OUTPUT
%   FFTPLOT     - Figure object; contains subplots for each input stream's
%                 FFT power spectrum.
%
% EXAMPLE USAGE:
%   % Plot FFT power spectra for session 1 with a frequency cutoff of 50 Hz
%   plotFFTpower(data, 1, 'Subject 427 - Session 1', 'fs', 'xmax', 50);
%
%   See also: plotFFTpower, plotstreamFFTmag, plotFFTmag
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

    %% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters_dev(mfilename); % For more details on default parameter values, see help configDefaultParameters.

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
    disp(['PLOTSTREAMFFTPOWER: Plotting FFT power spectrum for file: ',num2str(fileindex)])
    disp([' Stream: ',streamfieldname])
    
    if isempty(params.plotfilepath) & params.saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        error('SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

%% Prep FFTs
    [FFT,F] = preparestreamFFT(data(fileindex).(streamfieldname),data(fileindex).(fsfieldname));

    FFTPower = FFT.^2;

%% Prep colors
    FFTcolor = '#4CBB17';
    
%% Prep axis variables

    if strcmp('actual', params.xmax)
        params.xmax = max(F);
        disp(['  X axis (frequency) max set to actual max: ',num2str(params.xmax),' hz'])
    end

    currxlength = sum(F<params.xmax);
    currxticksize = floor(params.xmax/20); % Find total number of minutes per session - helper variable to determine ticks
    currxticks = 0:currxticksize:params.xmax;

    fftymax = 10^1;
    fftymin = 10^-20;
    fftyticks = [10^-20, 10^-10, 10^0];

%% Plot traces
    close all

    % Create tiled layout
    FFTplot = figure();

    % Plot stream FFT Power
    hold on;
    plot(F(1:currxlength), FFTPower(1:currxlength), 'Color', FFTcolor);
    set(gca, 'YScale', 'log')
    xlim([-0.1 params.xmax]);
    xticks(currxticks);
    ylim([fftymin fftymax]);
    yticks(fftyticks);
    title(streamfieldname);
    xlabel('Frequency (Hz)');
    ylabel('Power (Log)');
    set(gca,'TickDir','out')
    title(maintitle)
    hold off;


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
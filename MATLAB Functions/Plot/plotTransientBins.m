function [allbins] = plotTransientBins(data,whichfile,whichstream,whichtransients,maintitle,varargin)
% PLOTTRANSIENTBINS     Plots every bin for a session with detected
%                       transient peaks marked by circles.
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the
%                       streams to be plotted.
%
%       WHICHFILE:      The file number to plot (this can be set in a for
%                       loop to plot all files).
%
%       WHICHSTREAM:    String; The name of the field containing the stream
%                       used for transient detection. For example,
%                       'sigfiltz_normsession'.
%
%       WHICHTRANSIENTS: String; The name of the field containing the table
%                       of detection transients output by the function
%                       FINDSESSIONTRANSIENTS. For example,
%                       'sessiontransients_blmean_threshold3SD'.
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
%                       'C:\Users\rmdon\Box\Injection Transients\Figures\BinTraces_427_Morphine.png'
%
% OUTPUT:
%       ALLBINS:      A plot object containing subplots for each input
%                       stream.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Inputs
    inputs = struct(...
        'whichstream',[],...
        'whichtransients',[],...
        'saveoutput',[],...
        'plotfilepath',[]);
    inputs = parseArgsLite(varargin,inputs);

    % Prepare defaults and check for optional inputs
    inputs.whichstream = whichstream;
    inputs.whichtransients = whichtransients;

    if isempty(inputs.saveoutput)
        saveoutput = 0; % Defaults to 0 - skip saving plots
        inputs.saveoutput = saveoutput;
    else
        saveoutput = inputs.saveoutput;
    end


    disp(append('PLOTTRANSIENTBINS: Plotting bin traces with detected transients for file: ',num2str(whichfile)))
    disp(inputs)

    if isempty(inputs.plotfilepath) & saveoutput == 1 % If saveoutput is set to 1, plotfilepath is required
        disp('   ERROR: SAVEOUTPUT set to 1 but no PLOTFILEPATH specified. Provide PLOTFILEPATH or set SAVEOUTPUT to 0.')
    end

%% Prep colors
    tracecolor = '#4CBB17';
    transientcolor = '#0004ff';
    
%% Prep axis variables
    binsamples = data(whichfile).(whichtransients).BinSettings.binlengthsamples;
    binlengthmins = data(whichfile).(whichtransients).BinSettings.binlengthmins;
    nbins = data(whichfile).(whichtransients).BinSettings.nbins;
    
    binymax = ceil(max(data(whichfile).(whichstream))+(0.1*max(data(whichfile).(whichstream))));
    binymin = floor(min(data(whichfile).(whichstream))-(0.1*min(data(whichfile).(whichstream))));    
    binyticksize = round((binymax-binymin)/4,0); % Find size of ticks to generate 5 y axis ticks total
    binyticks = binymin:binyticksize:binymax;


%% Plot bin traces with transients
    close all

    % Create tiled layout
    allbins = tiledlayout(ceil(nbins/5), 5, 'Padding','compact', 'TileSpacing','compact');

     for bin = 0:(nbins-1)

        nexttile
        hold on
        
        % Prepare start and end indices for stream trace
        startbin = (bin*binsamples)+1;
        if length(data(whichfile).(whichstream)) > ((bin+1)*binsamples)
            endbin = (bin+1)*binsamples;
        else
            endbin = length(data(whichfile).(whichstream));
        end
        
        % Prepare transient max locations
        pklocs = data(whichfile).(whichtransients).transientquantification.maxloc(data(whichfile).(whichtransients).transientquantification.(append('Bin_',num2str(binlengthmins))) == (bin+1));

        % Plot stream trace with detected transients
        plot(data(whichfile).(whichstream)(startbin:endbin), 'Color', tracecolor)
        plot(pklocs-startbin,data(whichfile).(whichstream)(pklocs),'o', 'Color', transientcolor)
        axis([0 binsamples binymin binymax])
        xticks([(0:binlengthmins).*(binsamples/binlengthmins)])
        xticklabels((0:binlengthmins))
        yticks(binyticks)
        xlabel('Minute')
        title(append('Bin ',num2str(bin+1)))
        hold off
     end
    
    % Add a main title for the entire tiled layout
    title(allbins, maintitle, 'Interpreter', 'none');

    if saveoutput == 1
        set(gcf, 'Units', 'inches', 'Position', [0, 0, 8, 1.75*ntraces]);
        exportgraphics(gcf,append(plotfilepath, '.png'),'Resolution',300)
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
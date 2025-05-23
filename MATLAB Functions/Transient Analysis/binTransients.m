function [transientdata] = binTransients(transientdata,varargin)
% BINTRANSIENTS  Assigns each transient to a time bin within the session.
%
%   BINTRANSIENTS(DATA, STREAMFIELDNAME, FSFIELDNAME, TRANSIENTSFIELDNAME, 'PARAM1', VAL1, ...)
%   adds a 'Bin' variable to the transient quantification table, assigning
%   each transient to a bin based on its occurrence within the session.
%
% REQUIRED INPUTS:
%       TRANSIENTDATA       - Structure array output from findTransients.
%                             Must contain 'params.findTransients' and
%                             'transientquantification' fields.
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%       'binlengthmins'        - Numeric; length of each bin in minutes.
%                                Default: 5.
%
%       'nbinsoverride'        - Numeric; manual override to set the number of
%                                bins. If set to a value other than 0, this
%                                overrides the stream-length-based calculation
%                                of the number of bins per session.
%                                Default: 0.
%
%       'manuallydefinebins'   - Logical; Set to true (1) to input manually defined bin
%                                start and end indexes, rather than defining bins by time.
%                                Default: false (0).
%
%       'binstartfieldname'    - String; name of the field in DATA that contains
%                                custom bin start indexes. Only use if
%                                'manuallydefinebins' is set to true.
%
%       'binendfieldname'      - String; name of the field in DATA that contains
%                                custom bin end indexes. Only use if
%                                'manuallydefinebins' is set to true.
%
% OUTPUTS:
%   DATA                - Original transientdata structure with bins added to the
%                         specified table of transients.
%
% EXAMPLE:
%   transientdata = binTransients(transientdata, 'binlengthmins', 10);
%
% See also: findTransients
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
    addParameter(p, 'binlengthmins', defaultparameters.binlengthmins, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addParameter(p, 'nbinsoverride', defaultparameters.nbinsoverride, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
    addParameter(p, 'manuallydefinebins', defaultparameters.manuallydefinebins, @(x) islogical(x) || (isnumeric(x) && ismember(x, [0, 1])));
    addParameter(p, 'binstartfieldname', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'binendfieldname', '', @(x) ischar(x) || isstring(x));

    parse(p, varargin{:});

    % Retrieve all parsed inputs into params structure
    allparams = p.Results;

    % Initialize usedParams with always-used fields
    params = struct();
    params.binlengthmins = allparams.binlengthmins;
    params.nbinsoverride = allparams.nbinsoverride;
    params.manuallydefinebins = allparams.manuallydefinebins;
    
    % Conditionally add 'frequency' baqscalingtype specific params
    if params.manuallydefinebins==true
        params.binlengthmins = NaN;
        params.binstartfieldname = allparams.binstartfieldname;
        params.binendfieldname = allparams.binendfieldname;
    end
    
    if params.nbinsoverride > 0 && params.manuallydefinebins == true
        error('BIN DEFINEMENT ERROR - check optional inputs nbinsoverride and manuallydefinebins for conflicts. Note that only one of the two can be used at a time.')
    end
    
    % Display
    switch params.manuallydefinebins
        case false
            if params.nbinsoverride == 0
                disp(['BINTRANSIENTS: Add bin variable to transient quantification table. Binning transients into ', num2str(params.binlengthmins), ' minute bins.']) % Display bin length
                disp('   PARAMETERS:') % Display all parameters
                disp(params)
            else
                disp(['BINTRANSIENTS: Add bin variable to transient quantification table. Binning transients into ', num2str(params.nbinsoverride), ' total bins.']) % Display bin length
                disp('   PARAMETERS:') % Display all parameters
                disp(params)
            end
        case true
            disp(['BINTRANSIENTS: Add bin variable to transient quantification table. Binning transients by manually defined start indexes in ', params.binstartfieldname, ' and end indexes in ', params.binendfieldname, '.']) % Display bin length
            disp('   PARAMETERS:') % Display all parameters
            disp(params)
    end

    %% Bin Transients
    for eachfile = 1:length(transientdata)
        disp(['   Current File: ',num2str(eachfile)]) % Display which file is being processed

        fs = transientdata(eachfile).params.findTransients.fs;
        streamtotalsamples = transientdata(eachfile).params.findTransients.streamtotalsamples;
       
        binsamples = ceil(fs*60*params.binlengthmins); % Pull out n samples per bin

        if params.nbinsoverride == 0 && params.manuallydefinebins == false
            nbins = ceil(streamtotalsamples/binsamples); % Determine number of bins
            disp(['     ',num2str(nbins),' bins']) % Display which file is being processed
        elseif params.nbinsoverride > 0 && params.manuallydefinebins == false
            nbins = params.nbinsoverride;
            disp(['     NUMBER OF BINS MANUALLY SET: ', num2str(nbins),' bins']) % Display which file is being processed
        elseif params.manuallydefinebins == true
            binstartindexes = [transientdata(eachfile).(params.binstartfieldname)];
            binendindexes = transientdata(eachfile).(params.binendfieldname);
            nbins = length(binstartindexes);
            disp(['     BINS MANUALLY DEFINED: ', num2str(nbins),' bins']) % Display which file is being processed
        else
            error('BIN DEFINEMENT ERROR - check optional inputs nbinsoverride and manuallydefinebins for conflicts. Note that only one of the two can be used at a time.')
        end

        if params.manuallydefinebins == false
            transientdata(eachfile).transientquantification.(append('Bin_',num2str(params.binlengthmins),'min'))(1:height(transientdata(eachfile).transientquantification)) = NaN; % Clean the bin column to prevent errors if function parameters are changed and function is re-run
            for eachbin = 0:(nbins-1) % Add bin to transients column
                startbin = (eachbin*binsamples)+1; % Find bin start index
                endbin = (eachbin+1)*binsamples; % Find bin end index
                transientdata(eachfile).transientquantification.(append('Bin_',num2str(params.binlengthmins),'min'))((transientdata(eachfile).transientquantification.maxloc >= startbin & ...
                   transientdata(eachfile).transientquantification.maxloc < endbin)) = eachbin+1; % Add variable for bin to the transient quantification table. Transients are grouped by the bin start and stop indexes. Output column is labeled 'Bin_BINLENGTHMINS'
            end
            % Add bin transients settings to the transients table
            transientdata(eachfile).params.(mfilename).(append('Bin_',num2str(params.binlengthmins),'min')).nbins = nbins; % Add number of bins to settings table
            transientdata(eachfile).params.(mfilename).(append('Bin_',num2str(params.binlengthmins),'min')).binlengthmins = params.binlengthmins; % Add bin length in minutes to settings table        
            transientdata(eachfile).params.(mfilename).(append('Bin_',num2str(params.binlengthmins),'min')).binlengthsamples = binsamples; % Add bin length in samples to settings table
        else
            transientdata(eachfile).transientquantification.('Bin_Custom')(1:height(transientdata(eachfile).transientquantification)) = NaN; % Clean the bin column to prevent errors if function parameters are changed and function is re-run
            for eachbin = 1:nbins % Add bin to transients column
                startbin = binstartindexes(eachbin); % Find bin start index
                endbin = binendindexes(eachbin); % Find bin end index
                transientdata(eachfile).transientquantification.('Bin_Custom')((transientdata(eachfile).transientquantification.maxloc >= startbin & ...
                   transientdata(eachfile).transientquantification.maxloc < endbin)) = eachbin; % Add variable for bin to the transient quantification table. Transients are grouped by the bin start and stop indexes. Output column is labeled 'Bin_BINLENGTHMINS'
            end
            % Add bin transients settings to the transients table
            transientdata(eachfile).params.(mfilename).('Bin_Custom').nbins = nbins; % Add number of bins to settings table
            transientdata(eachfile).params.(mfilename).('Bin_Custom').binlengthmins = params.binlengthmins; % Add bin length in minutes to settings table        
            transientdata(eachfile).params.(mfilename).('Bin_Custom').binlengthsamples = binsamples; % Add bin length in samples to settings table
        end
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
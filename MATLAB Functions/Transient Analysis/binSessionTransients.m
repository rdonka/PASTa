function [data] = binSessionTransients(data,whichstream,whichfs,whichtransients,varargin)
% BINSESSIONTRANSIENTS  Assigns each transient to a time bin within the session.
%
%   BINSESSIONTRANSIENTS(DATA, WHICHSTREAM, WHICHFS, WHICHTRANSIENTS, 'PARAM1', VAL1, ...)
%   adds a 'Bin' variable to the transient quantification table, assigning
%   each transient to a bin based on its occurrence within the session.
%
% REQUIRED INPUTS:
%   DATA                - Structure array; contains at least the field
%                         with the data stream from which transients were
%                         detected, the sampling rate, and the fields
%                         containing the transient data with a column of
%                         transient max indexes.
%
%   WHICHSTREAM         - String; name of the field containing the data
%                         stream input to the findSessionTransients function
%                         used to identify transients. For example, 'sigz_normsession'.
%
%   WHICHFS             - String; name of the field containing the sampling
%                         rate of the streams. For example, 'fs'.
%
%   WHICHTRANSIENTS     - String; name of the parent field containing the
%                         table of transients to assign bins to. For example,
%                         'sessiontransients_blmin_3SD'.
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%   'whichtransientstable' - String; name of the field within WHICHTRANSIENTS
%                            that contains the quantification of individual
%                            transient events. Specify this if not using the
%                            format output from the FINDSESSIONTRANSIENTS functions.
%                            Default: 'transientquantification'.
%
%   'whichmaxlocs'         - String; name of the field containing the
%                            transient max locations (indexes) relative to
%                            the whole session. Specify this if not using the
%                            format output from the FINDSESSIONTRANSIENTS functions.
%                            Default: 'maxloc'.
%
%   'binlengthmins'        - Numeric; length of each bin in minutes.
%                            Default: 5.
%
%   'nbinsoverride'        - Numeric; manual override to set the number of
%                            bins. If set to a value other than 0, this
%                            overrides the stream-length-based calculation
%                            of the number of bins per session.
%                            Default: 0.
%
% OUTPUTS:
%   DATA                - Original data structure with bins added to the
%                         specified table of transients.
%
% EXAMPLE:
%   data = binSessionTransients(data, 'sigz_normsession', 'fs', 'sessiontransients_blmin_3SD', ...
%       'binlengthmins', 10);
%
% See also: findSessionTransients
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
    addParameter(p, 'whichtransientstable', defaultparameters.whichtransientstable, @(x) ischar(x) || isstring(x));
    addParameter(p, 'whichmaxlocs', defaultparameters.whichmaxlocs, @(x) ischar(x) || isstring(x));
    addParameter(p, 'binlengthmins', defaultparameters.binlengthmins, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addParameter(p, 'nbinsoverride', defaultparameters.nbinsoverride, @(x) isnumeric(x) && isscalar(x) && (x >= 0));

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;
    
    % Display
    disp(['BINSESSIONTRANSIENTS: Add bin variable to transient quantification table. Binning transients from ', whichtransients, ' into ', num2str(params.binlengthmins), ' minute bins.']) % Display bin length
    disp('   PARAMETERS:') % Display all parameters
    disp(params)

    %% Bin Transients
    for eachfile = 1:length(data)
        disp(['   Current File: ',num2str(eachfile)]) % Display which file is being processed

        binsamples = ceil(data(eachfile).(whichfs)*60*params.binlengthmins); % Pull out n samples per bin

        if params.nbinsoverride == 0
            nbins = ceil(length(data(eachfile).(whichstream))/binsamples); % Determine number of bins
            disp(append('     ',num2str(nbins),' bins')) % Display which file is being processed
        else
            nbins = params.nbinsoverride;
            disp(append('     NUMBER OF BINS MANUALLY SET: ', num2str(nbins),' bins')) % Display which file is being processed
        end

        % Add bin transients settings to the transients table
        data(eachfile).(whichtransients).BinSettings.nbins = nbins; % Add number of bins to settings table
        data(eachfile).(whichtransients).BinSettings.nbinsoverride = params.nbinsoverride; % Add bin number override to settings table
        data(eachfile).(whichtransients).BinSettings.binlengthmins = params.binlengthmins; % Add bin length in minutes to settings table        
        data(eachfile).(whichtransients).BinSettings.binlengthsamples = binsamples; % Add bin length in samples to settings table

        for eachbin = 0:(nbins-1) % Add bin to transients column
            startbin = (eachbin*binsamples)+1; % Find bin start index
            endbin = (eachbin+1)*binsamples; % Find bin end index
            data(eachfile).(whichtransients).(params.whichtransientstable).(append('Bin_',num2str(params.binlengthmins)))((data(eachfile).(whichtransients).(params.whichtransientstable).(params.whichmaxlocs) >= startbin & ...
               data(eachfile).(whichtransients).(params.whichtransientstable).(params.whichmaxlocs) < endbin)) = eachbin+1; % Add variable for bin to the transient quantification table. Transients are grouped by the bin start and stop indexes. Output column is labeled 'Bin_BINLENGTHMINS'
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
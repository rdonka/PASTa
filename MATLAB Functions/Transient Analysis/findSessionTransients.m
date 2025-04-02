function [data] = findSessionTransients(data,whichbltype,whichstream,whichthreshold,whichfs,varargin)
% FINDSESSIONTRANSIENTS  Detects and quantifies transients in a data stream
%                        for the entire session using specified baseline methods.
%
%   FINDSESSIONTRANSIENTS(DATA, WHICHBLTYPE, WHICHSTREAM, WHICHTHRESHOLD, WHICHFS, 'PARAM1', PARAM1, ...)
%   analyzes the specified data stream to detect transient events based on
%   the provided baseline type and threshold. The function supports various
%   baseline methods and allows customization through optional parameters.
%
% REQUIRED INPUTS:
%       DATA            - Structure array; each element corresponds to a session
%                         and must contain the fields specified by WHICHSTREAM,
%                         WHICHTHRESHOLD, and WHICHFS.
%
%       WHICHBLTYPE     - String; specifies the baseline type for transient detection.
%                         Options include:
%                           'blmean'   : Baseline is the mean of the pre-transient window.
%                           'blmin'    : Baseline is the minimum value within the pre-transient window.
%                           'localmin' : Baseline is the local minimum directly preceding the transient.
%
%       WHICHSTREAM     - String; name of the field in DATA containing the data stream
%                         to be analyzed (e.g., 'sigfiltz_normsession').
%
%       WHICHTHRESHOLD  - String; name of the field in DATA containing the numeric
%                         threshold values for transient detection (e.g., 'threshold_3SD').
%                         Thresholds should be precomputed and typically set to 2-3 standard deviations.
%
%       WHICHFS         - String; name of the field in DATA containing the sampling rate (fs)
%                         of the data stream.
%
%   OPTIONAL INPUT NAME-VALUE PAIRS:
%       'preminstartms'           - Numeric; start time (ms) of the pre-transient baseline window.
%                                   Default: 800 ms.
%
%       'preminendms'             - Numeric; end time (ms) of the pre-transient baseline window.
%                                   Default: 100 ms.
%
%       'posttransientms'         - Numeric; duration (ms) after the transient peak for analysis.
%                                   Default: 2000 ms.
%
%       'compoundtransientwindowms' - Numeric; window (ms) to search before and after each event
%                                     for compound transients. Default: 2000 ms.
%
%       'quantificationheight'    - Numeric; height (as a fraction of peak amplitude) at which to
%                                   characterize rise time, fall time, peak width, and area under
%                                   the curve (AUC). Must be between 0 and 1. Default: 0.5.
%
%       'outputtransientdata'     - Logical; if true (1), outputs cut data streams for each transient
%                                   event. If false (0), skips this output. Default: true (1).
%
%   OUTPUTS:
%       DATA            - Structure array; each element corresponds to a session and includes
%                         the following added fields:
%                           - sessiontransients_<WHICHBLTYPE>_<THRESHOLDLABEL>: A structure containing:
%                               - params: Structure of input parameters used for transient detection.
%                               - transientquantification: Table of quantified variables for each transient,
%                                 including amplitude, rise time, fall time, width, and AUC.
%                               - transientstreamlocs: Table of pre-transient baseline, transient peak,
%                                 rise, and fall locations for each transient.
%                               - transientstreamdata: Table of cut data streams from baseline start to
%                                 the end of the post-transient period for each transient event.
%
%   EXAMPLE:
%       % Detect transients using 'blmean' baseline method with default parameters
%       data = findSessionTransients(data, 'blmean', 'sigfiltz_normsession', 'threshold_3SD', 'fs');
%
%   See also: findSessionTransients_blmin, findSessionTransients_blmean, findSessionTransients_localmin
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

    % Import required and optional params into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details

    % Add optional name-value pair arguments with validation
    addParameter(p, 'preminstartms', defaultparameters.preminstartms, @isnumeric); % preminstartms: Numeric (ms); pre-transient to start the baseline period
    addParameter(p, 'preminendms', defaultparameters.preminendms, @isnumeric); % preminendms: Numeric (ms); ms pre-transient to start the baseline period
    addParameter(p, 'posttransientms', defaultparameters.posttransientms, @isnumeric); % posttransientms: Numeric (ms); ms post-transient to use as the oost-transient fall window
    addParameter(p, 'compoundtransientwindowms',  defaultparameters.compoundtransientwindowms, @isnumeric);  % Numeric (ms); Window size to search before and after each event for compound transients
    addParameter(p, 'quantificationheight', defaultparameters.quantificationheight, @(x) isnumeric(x) && x >= 0 && x <= 1); % quantificationheight: Numeric between 0 and 1; Height for quantification of transients (rise/fall/AUC)
    addParameter(p, 'outputtransientdata', defaultparameters.outputtransientdata, @islogical); % outputtransientdata: Logical; If set to 1 (true), data streams for individual transients will be added to the data structure

    parse(p, varargin{:});

    % Retrieve parsed params into params structure
    params = p.Results;

%% Call findSessionTransients subfunctions: findSessionTransients_blmin, findSessionTransients_blmean, or findSessionTransients_localmin
    if strcmp(whichbltype,'blmin') == true
        % Display settings
        disp('FINDSESSIONTRANSIENTS: Peak baseline determined by minimum value in the specified baseline window. WHICHBLTYPE set to blmin')
        disp(['     Transient data will be added to data structure as sessiontransients_blmin_',whichthreshold,'.'])
        disp('   PARAMETERS:') % Display all input values
        disp(params)
        
        % Find transients
        [data] = findSessionTransients_blmin(data,whichstream,whichthreshold,whichfs,params.preminstartms,params.preminendms,params.posttransientms,params.compoundtransientwindowms,params.quantificationheight,params.outputtransientdata);
    elseif strcmp(whichbltype,'blmean') == true
        % Display settings
        disp('FINDSESSIONTRANSIENTS: Peak baseline determined by mean value of the specified baseline window. WHICHBLTYPE set to blmean')
        disp(['     Transient data will be added to data structure as sessiontransients_blmean_',whichthreshold,'.'])
        disp('   PARAMETERS:') % Display all input values
        disp(params)
        % Find transients
        [data] = findSessionTransients_blmean(data,whichstream,whichthreshold,whichfs,params.preminstartms,params.preminendms,params.posttransientms,params.compoundtransientwindowms,params.quantificationheight,params.outputtransientdata);
    elseif strcmp(whichbltype,'localmin') == true
        % Display settings
        disp('FIND SESSION TRANSIENTS: Peak baseline determined by last local minimum in the specified baseline window. WHICHBLTYPE set to localmin')
        disp(['     Transient data will be added to data structure as sessiontransients_localmin_',whichthreshold,'.'])
        disp('   PARAMETERS:') % Display all input values
        disp(params)
        % Find transients
        [data] = findSessionTransients_localmin(data,whichstream,whichthreshold,whichfs,params.preminstartms,params.preminendms,params.posttransientms,params.compoundtransientwindowms,params.quantificationheight,params.outputtransientdata);
    else
        error('No viable transient baseline type specific. WHICHBLTYPE must be set to blmin, blmean, or localmin')
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
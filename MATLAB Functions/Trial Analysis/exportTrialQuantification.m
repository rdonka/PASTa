function [alltrialvals] = exportTrialQuantification(data,trialquantfieldname,exportfilepath,addvariablesfieldnames,varargin)
% EXPORTTRIALQUANTIFICATION  Compiles all trial quantificaiton values across sessions into a table and exports to a CSV file.
%
%   EXPORTTRIALQUANTIFICATION(DATA, TRIALQUANTFIELDNAME, EXPORTFILEPATH, ADDVARIABLESFIELDNAMES, 'PARAM1', VAL1, ...)
%   aggregates transient data from multiple sessions and exports the compiled table to a specified CSV file.
%
% REQUIRED INPUTS:
%   DATA                    - Struct array containing the output from QUANTIFYTRIALS 
%                             function with the name specified in TRIALQUANTFIELDNAME.
%
%   TRIALQUANTFIELDNAME     - String; name of the field produced by QUANTIFYTRIALS 
%                             function that contains the trial values to be exported.
%
%   EXPORTFILEPATH          - String; path to the folder where the CSV file will be
%                             saved. Note: The path must end with a forward slash.
%
%   ADDVARIABLESFIELDNAMES  - Cell array of strings; names of additional variables 
%                             from the data structure to include in the trial 
%                             quantification table. These variables will be added to 
%                             every row of the output table. At a minimum, this should 
%                             include the subject ID. If multiple sessions per subject 
%                             are included, ensure a session ID variable is also included.
%                               For example: {'Subject', 'SessionID', 'Treatment'}.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'exportfilename'    - String; custom name for the output CSV file. If not 
%                         specified, the file name will be  generated as 
%                         '<EXPORTFIELDNAME>_AllSessionExport_DAY-MONTH-YEAR.csv'.
%
% OUTPUTS:
%   ALLTRIALVALS   - Table; compiled table of all trial quantification values across all 
%                     sessions in the TRIALQUANTIFICATION field of the data structure. This table is also 
%                     saved as a CSV file at the specified export file path.
%
% EXAMPLE USAGE:
%   alltrialvals = exportTrialQuantification(data, trialquantfieldname, exportfieldname, exportfilepath, ...
%                       {'Subject', 'SessionID'}, 'exportfilename', 'transients_export.csv');
%
% Author:  Rachel Donka (2026)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Prepare Settings
    % Prepare default values
    defaultparameters = configDefaultParameters(mfilename); % For more details on default parameter values, see help configDefaultParameters.

    % Import required and optional inputs into a structure
    p = createParser(mfilename); % Create parser object with custom settings - see createParser helper function for more details
    addParameter(p, 'exportfilename', defaultparameters.exportfilename); % exportfilepath: string with the custom exportfilepath

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Set exportfilepath automatically if not specified
    if isempty(params.exportfilename)
        params.exportfilename = append(trialquantfieldname,'_AllTrialExport_',string(datetime("today")),'.csv');
    end
    
    % Display
    disp(['EXPORTTRIALQUANTIFICATION: Exporting all trial values as a csv file to: ', exportfilepath,'/',params.exportfilename]) % Display file path location
    disp(params)

    %% Prepare Trial Values for Export
    alltrialvals = table; % Prepare empty table

    for eachfile = 1:length(data)
        eachfiletrialvals = data(eachfile).(trialquantfieldname);
        
        for eachvariable = 1:length(addvariablesfieldnames)
            currvariable = char(addvariablesfieldnames(eachvariable));
            try 
                eachfiletrialvals.(currvariable)(1:height(eachfiletrialvals),1) = {data(eachfile).(currvariable)};
            catch
                disp(append('WARNING: File number ',num2str(eachfile), ' - failed to add variable: ', currvariable))
            end
        end

        if isempty(eachfiletrialvals)
            disp(append('   WARNING: Trial quantification table empty for file ', num2str(eachfile)))
        else
            alltrialvals = vertcat(alltrialvals,eachfiletrialvals);
        end
    end

    alltrialvals = movevars(alltrialvals,addvariablesfieldnames,"Before",1);
    writetable(alltrialvals,append(exportfilepath,params.exportfilename));
end

% Copyright (C) 2026 Rachel Donka
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
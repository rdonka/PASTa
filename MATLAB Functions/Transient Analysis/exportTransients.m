function [alltransients] = exportTransients(transientdata,exportfieldname,exportfilepath,addvariablesfieldnames,varargin)
% EXPORTTRANSIENTS  Compiles all transients across sessions into a table and exports to a CSV file.
%
%   EXPORTTRANSIENTS(TRANSIENTDATA, TRANSIENTQUANTIFICATIONFIELDNAME, EXPORTFILEPATH, ADDVARIABLESFIELDNAMES, 'PARAM1', VAL1, ...)
%   aggregates transient data from multiple sessions and exports the compiled table to a specified CSV file.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA     - Structure array of the output from FINDTRANSIENTS
%                       with the field you wish to export.
%
%   EXPORTFIELDNAME   - String; name of the field that contains the
%                     transients to be exported. 
%                     For individual transient events, use 'transientquantification'. 
%                     For session means of transient events, use 'transientsummary_session'. 
%                     For bin means of transient events, use 'transientsummary_<BINFIELDNAME>'.
%                     Ensure the TRANSIENTDATA structure contains the field
%                     you specify in EXPORTFIELDNAME.
%
%   EXPORTFILEPATH  - String; path to the folder where the CSV file will be
%                     saved. Note: The path must end with a forward slash.
%
%   ADDVARIABLESFIELDNAMES - Cell array of strings; names of additional variables 
%                     from the transientdata structure to include in the transients 
%                     table. These variables will be added to every row of 
%                     the output table. At a minimum, this should include 
%                     the subject ID. If multiple sessions per subject are 
%                     included, ensure a session ID variable is also included.
%                         For example: {'Subject', 'SessionID', 'Treatment'}.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'exportfilename'    - String; custom name for the output CSV file. If not 
%                         specified, the file name will be  generated as 
%                         '<EXPORTFIELDNAME>_AllSessionExport_DAY-MONTH-YEAR.csv'.
%
% OUTPUTS:
%   ALLTRANSIENTS   - Table; compiled table of all transients across all 
%                     sessions in the TRANSIENTDATA structure. This table is also 
%                     saved as a CSV file at the specified export file path.
%
% EXAMPLE USAGE:
%   alltransients = exportTransients(transientdata, exportfieldname, exportfilepath, ...
%                       {'Subject', 'SessionID'}, 'exportfilename', 'transients_export.csv');
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
    addParameter(p, 'exportfilename', defaultparameters.exportfilename); % exportfilepath: string with the custom exportfilepath

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Set exportfilepath automatically if not specified
    if isempty(params.exportfilename)
        params.exportfilename = append(exportfieldname,'_AllSessionExport_',string(datetime("today")),'.csv');
    end
    
    % Display
    disp(['EXPORTTRANSIENTS: Exporting transients as a csv file to: ', exportfilepath,'/',params.exportfilename]) % Display file path location
    disp(params)

    %% Prepare Transients for Export
    alltransients = table; % Prepare empty table

    for eachfile = 1:length(transientdata)

        eachfiletransients = transientdata(eachfile).(exportfieldname);
        
        for eachvariable = 1:length(addvariablesfieldnames)
            currvariable = char(addvariablesfieldnames(eachvariable));
            try 
                eachfiletransients.(currvariable)(1:height(eachfiletransients),1) = {transientdata(eachfile).(currvariable)};
            catch
                disp(append('WARNING: File number ',num2str(eachfile), ' - failed to add variable: ', currvariable))
            end
        end

        if isempty(eachfiletransients)
            disp(append('   WARNING: Transient quantification table empty for file ', num2str(eachfile)))
        else
            alltransients = vertcat(alltransients,eachfiletransients);
        end
    end

    alltransients = movevars(alltransients,addvariablesfieldnames,"Before",1);
    writetable(alltransients,append(exportfilepath,params.exportfilename));
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
function [alltransients] = exportTransients(data,whichtransients,exportfilepath,addvariables,varargin)
% EXPORTTRANSIENTS  Compiles all transients across sessions into a table and exports to a CSV file.
%
%   EXPORTTRANSIENTS(DATA, WHICHTRANSIENTS, EXPORTFILEPATH, ADDVARIABLES, 'PARAM1', VAL1, ...)
%   aggregates transient data from multiple sessions and exports the compiled table to a specified CSV file.
%
% REQUIRED INPUTS:
%   DATA            - Structure array; must contain at least the output 
%                     from FINDTRANSIENTS.
%
%   WHICHTRANSIENTS - String; name of the field containing the table of 
%                     transients to export (e.g., 'sessiontransients_blmin_3SD').
%
%   EXPORTFILEPATH  - String; path to the folder where the CSV file will be
%                     saved. Note: The path must end with a forward slash.
%
%   ADDVARIABLES    - Cell array of strings; names of additional variables 
%                     from the data structure to include in the transients 
%                     table. These variables will be added to every row of 
%                     the output table. At a minimum, this should include 
%                     the subject ID. If multiple sessions per subject are 
%                     included, ensure a session ID variable is also included.
%                         For example: {'Subject', 'SessionID', 'Treatment'}.
%
% OPTIONAL INPUT NAME-VALUE PAIR ARGUMENTS:
%   'whichtransientstable' - String; name of the field within WHICHTRANSIENTS 
%                            that contains the quantification of individual 
%                            transient events. This input only needs to be 
%                            specified if not using the format output from 
%                            the FINDSESSIONTRANSIENTS functions.
%                            Default: 'transientquantification'.
%
%   'filename'             - String; custom name for the output CSV file. 
%                            If not specified, the file name will be 
%                            generated as 'WHICHTRANSIENTS_AllSessionExport_DAY-MONTH-YEAR.csv'.
%
% OUTPUTS:
%   ALLTRANSIENTS   - Table; compiled table of all transients across all 
%                     sessions in the data structure. This table is also 
%                     saved as a CSV file at the specified export file path.
%
% EXAMPLE USAGE:
%   alltransients = exportTransients(data, 'sessiontransients_blmin_3SD', exportfilepath, {'Subject', 'SessionID'}, 'filename', 'transients_export.csv');
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
    addParameter(p, 'whichtransientstable', defaultparameters.whichtransientstable, @(x) ischar(x) || isstring(x)); % whichtransientstable: string with the name of the transients table
    addParameter(p, 'filename', defaultparameters.filename); % filename: string with the custom filename

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Set filename automatically if not specified
    if isempty(params.filename)
        params.filename = append(whichtransients,'_AllSessionExport_',string(datetime("today")),'.csv');
    end
    
    % Display
    disp(['EXPORTTRANSIENTS: Exporting transients from ', whichtransients, ' to a csv file.']) % Function display
    disp(['     File will be output to the folder location: ', exportfilepath]) % Display file path location

    disp('   PARAMETERS:') % Display all input values
    disp(params)

    %% Prepare Transients for Export
    alltransients = table; % Prepare empty table

    for eachfile = 1:length(data)

        eachfiletransients = data(eachfile).(whichtransients).(params.whichtransientstable);
        
        for eachvariable = 1:length(addvariables)
            currvariable = char(addvariables(eachvariable));
            try 
                eachfiletransients.(currvariable)(1:height(eachfiletransients),1) = {data(eachfile).(currvariable)};
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

    alltransients = movevars(alltransients,addvariables,"Before",1);
    writetable(alltransients,append(exportfilepath,params.filename));
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
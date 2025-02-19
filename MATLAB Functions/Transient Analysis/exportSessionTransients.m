function [alltransients] = exportSessionTransients(data,whichtransients,exportfilepath,addvariables,varargin)
% EXPORTSESSIONTRANSIENTS   Creates a table of all transients across all
%                           sessions in the data structure and exports to a
%                           csv file. Used in conjunction with
%                           FINDSESSIONTRANSIENTS.
%
% Copyright (C) 2025 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:               Data structure; This is a structure that contains 
%                           at least the output from FINDSESSIONTRANSIENTS
%
%       WHICHTRANSIENTS:    String; The name of the parent field containing 
%                           the table of transients that you want to export. 
%                           For example, 'sessiontransients_blmin_3SD'.
%
%       EXPORTFILEPATH:     String; Path to the folder location where the
%                           created table should be saved to. NOTE: The
%                           path must end in a forward slash.
%
%       ADDVARIABLES:       Cell array containing any additional variables
%                           from the data structure to be added to the
%                           transients table. Variables will be added to
%                           every row of the output structure. Cell array
%                           inputs must be the names of fields in the data
%                           structure. At a minimum, this should contain
%                           the subject ID. If multiple sessions per
%                           subject are included in the data structure,
%                           make sure a session ID variable is also
%                           included. 
%                           For example: {'Subject', 'SessionID', 'Treatment'}
%                       
% OPTIONAL INPUTS:
%       WHICHTRANSIENTSTABLE: String; The name of the field within WHICHTRANSIENTS
%                           that contains the quantification of individual 
%                           transient events. This input only needs to be
%                           specified if not using the format output from
%                           the FINDSESSIONTRANSIENTS functions.
%                           Default: 'transientquantification'.
%
%
%       FILENAME:           String; Custom name for the output csv file. If
%                           not specified, the file name will be generated
%                           as 'WHICHTRANSIENTS_AllSessionExport_DAY-MONTH-YEAR.csv'
%
% OUTPUTS:
%       This function outputs a csv file with all transients for all
%       sessions in the data structure. The file will be output at the
%       specified file path. If the function is called into an object, the
%       table ALLTRANSIENTS will also be saved to an object in the MATLAB
%       workspace. 
%       For example: alltransients = exportSessionTransients(...)
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichtransients',[],...
        'exportfilepath',[],...
        'addvariables',[],...
        'whichtransientstable', [],...
        'filename',[]);
    inputs = parseArgsLite(varargin,inputs);
    
    % Prepare defaults and check for optional inputs
    inputs.whichtransients = whichtransients;
    inputs.exportfilepath = exportfilepath;
    inputs.addvariables = addvariables;


    if isempty(inputs.whichtransientstable) % Field containing transient quantification
        whichtransientstable = 'transientquantification';
        inputs.whichtransientstable = whichtransientstable;
    else
        whichtransientstable = inputs.whichtransientstable;
    end

    if isempty(inputs.filename) % Field containing transient max location
        filename = append(whichtransients,'_AllSessionExport_',string(datetime("today")),'.csv');
        inputs.filename = filename;
    else
        filename = inputs.filename;
    end

    disp(append('EXPORT SESSION TRANSIENTS: Exporting transients from ', whichtransients, ' to a csv file.')) % Function display
    disp(append('     File will be output to the folder location: ', exportfilepath)) % Display file path location

    disp('INPUTS:') % Display all input values
    disp(inputs)

    %% Prepare Transients for Export
    alltransients = table; % Prepare empty table

    for eachfile = 1:length(data)

        eachfiletransients = data(eachfile).(whichtransients).(whichtransientstable);
        
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
    writetable(alltransients,append(exportfilepath,filename));
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
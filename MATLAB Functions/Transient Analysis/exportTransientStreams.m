function [alltransientstruct] = exportTransientStreams(transientdata,transientstreamfieldnames,addvariablesfieldnames,varargin)
% EXPORTTRANSIENTSTREAMS  Compiles all transients across sessions into a
% data structure. Primarily used for additional quantification and plotting
% of individual transient events.
%
%   EXPORTTRANSIENTSTREAMS(TRANSIENTDATA, TRANSIENTSTREAMFIELDNAMES, ADDVARIABLESFIELDNAMES)
%   aggregates transient data from multiple sessions. Each struct row
%   corresponds to one transient event.
%
% REQUIRED INPUTS:
%   TRANSIENTDATA               - Structure array of the output from FINDTRANSIENTS
%                                 containing the fields 'transientquantification' and
%                                 'transientstreamdata' at a minimum.
%
%   TRANSIENTSTREAMFIELDNAMES   - Cell array of strings; names of the fields    
%                                 that contain the cut transient data streams
%                                 to be exported.
%
%   ADDVARIABLESFIELDNAMES       - Cell array of strings; names of additional variables 
%                                  from the transientdata structure to include in the transients 
%                                  table. These variables will be added to every row of 
%                                  the output table. At a minimum, this should include 
%                                  the subject ID. If multiple sessions per subject are 
%                                  included, ensure a session ID variable is also included.
%                                  For example: {'Subject', 'SessionID', 'Treatment'}.
%
% OUTPUTS:
%   ALLTRANSIENTSSTRUCT   - Struct; compiled struct of all transients across all sessions in the 
%                           TRANSIENTDATA structure with all fields in 'transientquantification', 
%                           all fields in ADDVARIABLESFIELDNAMES, and each stream specified in 
%                           TRANSIENTSTREAMFIELDNAMES.
%
% EXAMPLE USAGE:
%   alltransientstruct = exportTransientStreams(transientdata, transientstreamfieldnames, addvariablesfieldnames);
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Display
    disp(['EXPORTTRANSIENTSTREAMS: Compiling structure of all transient event quantification and streams.'])

    %% Create struct with one row per transient
    alltransientstruct = struct([]); % Prepare empty struct

    for eachfile = 1:length(transientdata)
        disp(['Preparing File: ', num2str(eachfile)]);
        % Temp structure for current file
        eachfilestruct = struct([]);

        % Prepare field names for transient quantification table
        transientquantificationfieldnames = fieldnames(transientdata(eachfile).transientquantification);
        transientquantificationfieldnames = setdiff(transientquantificationfieldnames, {'Properties','Row','Variables'});
        
        for eachevent = 1:height(transientdata(eachfile).transientquantification)
            eacheventstruct = struct();
            % Add variables from ADDVARIABLESFIELDNAMES
            for eachvariable = 1:length(addvariablesfieldnames)
                currvariable = char(addvariablesfieldnames(eachvariable));
                try 
                    eacheventstruct.(currvariable) = transientdata(eachfile).(currvariable);
                catch
                    disp(append('WARNING: File number ',num2str(eachfile), ' - failed to add variable: ', currvariable))
                end
            end

            % Add variables from transient quantification table
            for eachtransientquantificationfieldname = 1:length(transientquantificationfieldnames)
                currtransientquantificationfield = transientquantificationfieldnames{eachtransientquantificationfieldname};
                eacheventstruct.(currtransientquantificationfield) = transientdata(eachfile).transientquantification.(currtransientquantificationfield)(eachevent);
            end
            
            % Add fields for each stream in TRANSIENTSTREAMFIELDNAMES
            for eachstreamfieldname = 1:length(transientstreamfieldnames)
                currstreamfieldname = transientstreamfieldnames{eachstreamfieldname};
                eacheventstruct.(currstreamfieldname) = transientdata(eachfile).(currstreamfieldname)(eachevent,:);
            end
            
            % Add event to temp current file structure
            if isempty(eachfilestruct)
                eachfilestruct = eacheventstruct; % seed with first
            else
                eachfilestruct(end+1,1) = eacheventstruct; % grow struct array
            end
        end

        % Combine file structures into overall structure output
        if isempty(alltransientstruct)
            alltransientstruct = eachfilestruct; % seed with first
        else
            alltransientstruct = [alltransientstruct; eachfilestruct];
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
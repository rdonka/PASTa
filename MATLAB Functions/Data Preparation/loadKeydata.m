function [data] = loadKeydata(experimentkey)
% LOADKEYDATA       Loads individual sessions into a main data structure.
%                   Individual sessions must have been previously extracted
%                   and saved as MATLAB data structures with one of the
%                   following functions: EXTRACTTDTDATA, EXTRACTDORICDATA,
%                   EXTRACTNPDATA, EXTRACTGENERICDATA
%
% Copyright (C) 2024 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       EXPERIMENTKEY:  Data Structure; A prepared data structure with at 
%                       minimum the field 'ExtractedFolderPath' containing
%                       the full path to the individual session data
%                       structures to be loaded. The EXPERIMENTKEY can be
%                       prepared with the LOADKEYS function.
%
% OUTPUTS:
%       DATA:           Data Structure; A data structure with each individual
%                       extracted session block as a row.
%
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

data = experimentkey; % Prepare data structure.

disp('LOADING DATA: Data must be previously extracted and saved as MATLAB data structures.')
    for eachfile = 1:length(experimentkey)
        try
            fprintf('Loading file number: %.f \n',eachfile) % Display which file is loading
            
            blockdata = load(strcat(data(eachfile).ExtractedFolderPath,'_extracted.mat')); % Load in pre-extracted block data structure
            
            blockdatafields = fieldnames(blockdata); % Find all the fieldnames of the pre-extracted block
    
            for eachfield = 1:length(blockdatafields) % Add each field to the main data structure
                data(eachfile).(blockdatafields{eachfield}) = blockdata.(blockdatafields{eachfield});
            end
        catch ME
            fprintf('ERROR: File %s\n', eachfile);
            fprintf('Message: %s\n', ME.message);
            continue;
        end
        fprintf('SUCCESS: File %.f \n',eachfile)
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
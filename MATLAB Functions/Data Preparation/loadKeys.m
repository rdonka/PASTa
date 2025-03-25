function [experimentkey] = loadKeys(rootdirectory,subjectkeyname,filekeyname)
% LOADKEYS  Combine subject key and file key into a single data structure, 
%           appending the provided computer user path to relevant fields.
%
%   [EXPERIMENTKEY] = LOADKEYS(ROOTDIRECTORY, SUBJECTKEYNAME, FILEKEYNAME)
%   reads two CSV files: a subject key and a file key, merges them on a
%   shared column (typically "SubjectID"), and then appends ROOTDIRECTORY
%   to the folder paths in the file key. This produces a single struct 
%   (EXPERIMENTKEY) ready for downstream analysis.
%
%   REQUIRED INPUTS:
%       ROOTDIRECTORY     - String. The top-level directory path unique to
%                           your system. For example: 'C:\Users\rmdon\'.
%                           Must end with a slash/backslash.
%
%       SUBJECTKEYNAME    - String. The name (or full path) of the subject
%                           key CSV file (e.g. 'mySubjectKey.csv'). Must
%                           contain at least the variable 'Subject'.
%                           If empty (''), the subject key is skipped and
%                           only FILEKEYNAME is loaded.
%
%       FILEKEYNAME       - String. The name (or full path) of the file key
%                           CSV file (e.g. 'myFileKey.csv'). Must contain 
%                           'Subject', 'Folder', 'RawFolderPath', and 
%                           'ExtractedFolderPath' columns at minimum.
%                           Paths in 'RawFolderPath' and 'ExtractedFolderPath'
%                           should each end with a slash/backslash.
%
% OUTPUTS:
%       EXPERIMENTKEY     - A struct array combining file and subject key
%                           data. The 'rootdirectory' is prepended to 
%                           each row's 'RawFolderPath' and 'ExtractedFolderPath',
%                           while the 'Folder' name is appended to both. 
%
%   EXAMPLE:
%       % Suppose your top-level path is: rootdirectory = 'C:\Users\rmdon\';
%       % You have the subject key and file key CSVs:
%       subjKey = 'subjectKey.csv';
%       fileKey = 'fileKey.csv';
%
%       % Load them into one experiment key structure:
%       experimentkey = loadKeys(compPath, subjKey, fileKey);
%
% Author:  Rachel Donka (2025)
% License: GNU General Public License v3. See end of file for details.
% Stored in the PASTa GitHub Repository: https://github.com/rdonka/PASTa
% For detailed instructions, see the PASTa user guide: https://rdonka.github.io/PASTaUserGuide/

%% Load the Subject Key and File Key CSV files as tables
% If 'subjectkeyname' is non-empty, merge the subject and file keys. Otherwise, just load the file key.

    if strcmp(subjectkeyname, "") == 0 % Subject key is provided, so read both and then attempt to join
        subjectkey = readtable(subjectkeyname, 'Decimal',',', 'Delimiter',','); % Load subject key
        filekey = readtable(filekeyname, 'Decimal',',', 'Delimiter',','); % Load file key
        try 
            % Merge the file key and subject key into one table, then convert to struct
            [experimentkey] = table2struct(join(filekey, subjectkey));
        catch
            % If join fails, likely the subject IDs in file key don't appear in subject key
            fprintf('ERROR: Could not join file key [%s] and subject key [%s].\n', filekeyname, subjectkeyname);
            fprintf('       Check that all "Subject" values in the file key exist in the subject key.\n\n');
            disp('Unique Subject IDs in Subject Key:');
            disp(unique(subjectkey.SubjectID));
            disp('Unique Subject IDs in File Key:');
            disp(unique(filekey.SubjectID));

            % Return error
            error('loadKeys:JoinError', 'Failed to join the subject and file key tables.');
        end
    else % No subject key provided, so only read the file key
        filekey = readtable(filekeyname, 'Decimal',',', 'Delimiter',',');
        [experimentkey] = table2struct(filekey);
    end
    
    % Append the root directory and folder name to each record.
    % The 'Folder' field is appended to the end of 'RawFolderPath' and 'ExtractedFolderPath', while ROOTDIRECTORY is prepended to both
    for eachfile = 1:length(experimentkey) % Prepend user path and append folder name
        [experimentkey(eachfile).RawFolderPath] = strcat(rootdirectory, experimentkey(eachfile).RawFolderPath, experimentkey(eachfile).Folder);
        [experimentkey(eachfile).ExtractedFolderPath] = strcat(rootdirectory, experimentkey(eachfile).ExtractedFolderPath, experimentkey(eachfile).Folder);
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

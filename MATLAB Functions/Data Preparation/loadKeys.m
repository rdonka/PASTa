function [experimentkey] = loadKeys(computeruserpath,subjectkeyname,filekeyname)
% LOADKEYS    Combines subject key and file key into a data structure, and 
%             appends the provided computeruserpath to the paths in the
%             file key.
%
% INPUTS:
%       COMPUTERUSERPATH:   String; A variable containing the unique portion
%                           of the path to the data folders for the users
%                           specific computer. For example, 'C:\Users\rmdon\'. 
%                           NOTE: Make sure it ends in a forward slash.
%
%       SUBJECTKEYNAME:     String; A variable containing a string with the
%                           name of the subject key csv file for the experiment.
%                           Must contain the field 'Subject' at a minumum.
%                           NOTE: If subjectkeyname is left empty (set to
%                           ''), the subject key will be skipped and only 
%                           the file key will be loaded into the experiment 
%                           key output.
%
%       FILEKEYNAME:        String; A variable containing a string with the
%                           name of the file key csv file for the experiment.
%                           Must contain the fields 'Subject', 'Folder',
%                           'RawFolderPath', and 'ExtractedFolderPath', at
%                           a minimum. 
%                           NOTE: Folder paths in 'RawFolderPath' and
%                           'ExtractedFolderPath' must end with a '\'.
%
% OUTPUTS:
%       EXPERIMENTKEY:      A data structure that includes the joined file 
%                           key and subject key with computer user path 
%                           appended to the front and the folder name 
%                           appended to the end of the RawFolderPath and 
%                           ExtractedFolderPath.
%
% Written by R M Donka, August 2023
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

    %% Load in the subject key and file key csv files as tables
    if strcmp(subjectkeyname, "") == 0 % If the subject key name is not empty, match it to the file key 
        subjectkey = readtable(subjectkeyname, 'Decimal',',', 'Delimiter',','); % Load subject key
        filekey = readtable(filekeyname, 'Decimal',',', 'Delimiter',','); % Load file key
        try
            [experimentkey] = table2struct(join(filekey, subjectkey)); % Merge the subject and file keys into data structure called experimentkey
        catch
            fprintf(1,'ERROR joining file key and subject key. The key variable for the right table must contain all values in the key variable for the left table.');
            disp(' ');
            disp('Check that all subject IDs contained in the file key are contained in the subject key.');
            disp('Unique Subject IDs in Subject Key:');
            disp(unique(subjectkey.Subject));
            disp('Unique Subject IDs in File Key:');
            disp(unique(filekey.Subject));
        end
    else % If the subejct key name is empty, only load the file key    
        filekey = readtable(filekeyname, 'Decimal',',', 'Delimiter',',');
        [experimentkey] = table2struct(filekey);
    end
    
    % Add full filepath to the experiment key by combining the parent folder path and the specific computer user path
    for eachfile = 1:length(experimentkey)
        [experimentkey(eachfile).RawFolderPath] = strcat(computeruserpath, experimentkey(eachfile).RawFolderPath, experimentkey(eachfile).Folder); % Raw file path
        [experimentkey(eachfile).ExtractedFolderPath] = strcat(computeruserpath, experimentkey(eachfile).ExtractedFolderPath, experimentkey(eachfile).Folder); % Output path for the figures
    end
end
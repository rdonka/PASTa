function [data] = loadKeydata(experimentkey)
% LOADFPDATA    Loads in pre-extracted subject dataframes and combines it into one dataframe.
% Requires data to be pre-extracted and saved into structures by EXTRACTFPDATA.
%
% INPUTS:
%       EXPERIMENTKEY:  A prepared data structure with at minimum the 
%                       ExtractedFolderPath to locate individual block
%                       structures.
%
% OUTPUTS:
%       DATA:           A data structure called "data" with each individual
%                       extracted block as a row.
%
% Written by R M Donka, August 2023
% Stored in RoitmanPhotometry GitHub repository, see Wiki for additional notes.

data = experimentkey; % Prepare data structure.

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


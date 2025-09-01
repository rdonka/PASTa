function [] = extractTDTdata_multiblock(experimentkey,sigstreamnamefield,baqstreamnamefield,varargin)
% EXTRACTTDTDATA  Extract and save fiber photometry data from TDT blocks
% for sessions where multiple subjects are saved to the same block.
%                 (Synapse, Tucker Davis Technologies).
%
%   EXTRACTTDTDATA(EXPERIMENTKEY, SIGSTREAMFIELDNAME, BAQSTREAMFIELDNAME, 
%   'PARAM1', VAL1, 'PARAM2', VAL2, ...) reads raw TDT block 
%   data from the folders specified in the the 'RawFolderPath' field in the 
%   experiment key, extracts signal and background streams from the stored
%   listing stream fields specified by 'sigstreamname' and 'baqstreamname', 
%   trims the data, and saves the results as MATLAB .mat files to the path 
%   specified in the 'ExtractedFolderPath' field in the experiment key with 
%   the SubjectID appended in the file name.
%
%   For each TDT block, the matching "signal" stream is selected from 
%   the experiment key field matching the SIGSTREAMFIELDNAME, and the 
%   matching "background" stream is selected from the experiment key field 
%   matching the BAQSTREAMFIELDNAME. Field names must be provided and
%   matched to the proper SubjectID in the file key.
%   The final extracted data structure is saved as: 
%         <ExtractedFolderPath>_<SubjectID>_extracted.mat
%
% REQUIRED INPUTS:
%       EXPERIMENTKEY:  Data Structure; A prepared data structure with at 
%                       minimum the fields:
%                         'RawFolderPath': field containing the full path to 
%                           the raw data block containing the session data.
%                         'ExtractedFolderPath': field containing the full
%                           path to the location where the extracted data
%                           block containing the session data should be
%                           output.
%                         'sigstreamfield': field containing the name of
%                           the stream containing the signal field for the
%                           session (e.g., 'x465A' or 'x465B')
%                         'baqstreamfield': field containing the name of
%                           the stream containing the background field for
%                           the session (e.g., 'x405A' or 'x405B')
%                       The EXPERIMENTKEY can be prepared with the LOADKEYS function.
%
%       SIGSTREAMFIELDNAME:  - String; Field name for the stream name of 
%                               the signal field of each session.
%                               Example: 'StoresListing_465'
%
%       BAQSTREAMFIELDNAME:  - String; Field name for the stream name of 
%                               the background field of each session.
%                               Example: 'StoresListing_465'
%
%  OPTIONAL INPUT NAME-VALUE PAIRS:
%       'trim'        - Numeric; Number of seconds to trim from the start 
%                       and end of each recording. Default: 5
% 
%       'skipexisting'-  Numeric (0 or 1); If 1, skip extracting any
%                       session for which an output file already exists.
%                       If 0, re-extract and overwrite. Default: 1
%
%   EXAMPLE:
%       extractTDTdata_multiblock(experimentkey, 'StoresListing_465', ...
%           'StoresListing_405', 'trim', 10, 'skipexisting', 0);
%
%   See also: TDTbin2mat
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
    addParameter(p, 'trim', defaultparameters.trim, @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'})); % trim: input must be a numeric value, nonnegative, and integer
    addParameter(p, 'skipexisting', defaultparameters.skipexisting, @(x) any(x == [0 1])); % skipexisting: input must be 0 or 1
    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Main display and function inputs
    if params.skipexisting == 0
        disp('EXTRACTTDTDATA: all blocks will be extracted.')
    elseif params.skipexisting == 1
        disp('EXTRACTTDTDATA: pre-extracted blocks will be skipped.')
    end

    disp('  PARAMETERS:') % Display all input parameters values
    disp(params)

    disp('    Parameters using default values:'); % Display input parameters values set to default values
    disp(p.UsingDefaults);

    %% Extract Raw Data
    for eachfile = 1:length(experimentkey)
        % Prepare a new block data structure for this file
        blockdata = struct();
        if isnumeric(experimentkey(eachfile).SubjectID) % Identify SubjectID and convert to string if needed
            blockdata.SubjectID = num2str(experimentkey(eachfile).SubjectID);
        else
            blockdata.SubjectID = experimentkey(eachfile).SubjectID;
        end
        blockdata.RawFolderPath = char(experimentkey(eachfile).RawFolderPath);
        blockdata.ExtractedFolderPath = char(experimentkey(eachfile).ExtractedFolderPath);

        % Construct output file name (where the extracted data will be saved)
        extractedMatFile = strcat(blockdata.ExtractedFolderPath,'_',blockdata.SubjectID,'_extracted.mat');
        
        % If set to skip existing files, check whether the output .mat already exists
        if params.skipexisting == 1 && isfile(extractedMatFile)
            fprintf('Skipping file #%d (already exists): %s\n', eachfile, extractedMatFile);
            continue
        else
            % Attempt to load and process the TDT block
            try
                fprintf('Extracting file #%d: %s\n', eachfile, blockdata.RawFolderPath); % Display which file is loading

                % Read in the entire TDT block data using TDTbin2mat - see function for details.
                alldata=TDTbin2mat(blockdata.RawFolderPath);
        
                % Identify all stream and epoc fields
                streams = string(fieldnames(alldata.streams)); % Find all stream names for the current block
                epocs = string(fieldnames(alldata.epocs)); % Find all epoch names for the current block
                
                % Find the signal and background stream names for the current session
                currsigname = experimentkey(eachfile).(sigstreamnamefield); % Extract the name of the signal field for the current file
                currbaqname = experimentkey(eachfile).(baqstreamnamefield); % Extract the name of the background field for the current file

                % Display warning if no matching signal or background was found
                if ~isfield(alldata.streams,currsigname)
                    warning('No valid SIGNAL stream found for file #%d.', eachfile);
                end
                if ~isfield(alldata.streams,currbaqname)
                    warning('No valid BACKGROUND stream found for file #%d.', eachfile);
                end
                             
                disp(append('     Signal stream: ',currsigname)) % Display current file signal field name
                disp(append('     Background stream: ',currbaqname)) % Display current file background field name

                % Pull out timing/metadata
                blockdata.date = alldata.info.date; % Add date to data 
                blockdata.sessionduration = alldata.info.duration; % Add session duration to data
                blockdata.starttime = alldata.info.utcStartTime;
                blockdata.stoptime = alldata.info.utcStopTime;
                
                % Add the sampling rate for the signal stream
                blockdata.fs = alldata.streams.(currsigname).fs;

                % Trim the data at the beginning and end, as specified
                trimsamples = round(params.trim*blockdata.fs); % Find number of samples to trim from beginning and end of session
                blockdata.sig = double(alldata.streams.(currsigname).data((trimsamples+1):end-trimsamples)); % Add signal to the data structure and trim
                blockdata.baq = double(alldata.streams.(currbaqname).data((trimsamples+1):end-trimsamples)); % Add background to the data structure and trim

                % Adjust epoch event times for the trimmed data based on epoc onset
                for eachepoc = 1:length(epocs)
                    thisepoc = char(epocs(eachepoc));
                    blockdata.(thisepoc) = round(((alldata.epocs.(thisepoc).onset-alldata.time_ranges(1)).*blockdata.fs)-(trimsamples));
                end

                % Add a record of which function + parameters were used
                blockdata.params.(mfilename) = params;
                blockdata.params.(mfilename).extractFunction = mfilename;
            
                % Save the extracted data struct as <ExtractedFolderPath>_extracted.mat
                save(extractedMatFile,'-struct', 'blockdata');

            catch ME % % If there's an error (e.g., path not found), log it and continue.
                fprintf('ERROR extracting file #%d. Verify raw and extract file paths.\n', eachfile);
                fprintf('   Raw path: %s\n', blockdata.RawFolderPath);
                fprintf('   Sig field name: %s\n', currsigname);
                fprintf('   Baq field name: %s\n', currbaqname);
                fprintf('   Error ID: %s\n', ME.identifier);
                fprintf('   Message: %s\n', ME.message);
                continue
            end
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
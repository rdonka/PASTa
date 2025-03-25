function [] = extractTDTdata(rawfolderpaths,extractedfolderpaths,sigstreamnames,baqstreamnames,varargin)
% EXTRACTTDTDATA  Extract and save fiber photometry data from TDT blocks
%                 (Synapse, Tucker Davis Technologies).
%
%   EXTRACTTDTDATA(RAWFOLDERPATHS, EXTRACTEDFOLDERPATHS, SIGSTREAMNAMES, 
%   BAQSTREAMNAMES, 'PARAM1', VAL1, 'PARAM2', VAL2, ...) reads raw TDT block 
%   data from the folders specified in RAWFOLDERPATHS, extracts signal and 
%   background streams, trims the data, and saves the results to 
%   EXTRACTEDFOLDERPATHS in MATLAB .mat files.
%
%   For each TDT block, the matching "signal" stream is selected from 
%   SIGSTREAMNAMES, and the matching "background" stream is selected from 
%   BAQSTREAMNAMES. If multiple possible stream names exist across session,
%   all possible names should be provided in SIGSTREAMNAMES and 
%   BAQSTREAMNAMES. The final extracted data structure is saved as: 
%         <ExtractedFolderPath>_extracted.mat
%
% REQUIRED INPUTS:
%       RAWFOLDERPATHS        - String array of raw TDT block folder paths.
%                               The string array should contain one column 
%                               with each full path in a separate row. If 
%                               using the loadKeys function, this can be 
%                               created from the experiment key. For example:
%                                   rawfolderpaths = string({experimentkey.RawFolderPath})';
%
%       EXTRACTEDFOLDERPATHS  - String array of corresponding output paths 
%                               where extracted data is saved. If using the 
%                               loadKeys function, this can be created from 
%                               the experiment key. For example: 
%                                   extractedfolderpaths = string({experimentkey.ExtractedFolderPath})';
%
%       SIGSTREAMNAMES        - Cell array of possible stream names for 
%                               the "signal" channel (e.g. {'x65A','465A'}). 
%                               NOTE: Only one stream per file can be
%                               treated as signal.
%
%       BAQSTREAMNAMES        - Cell array of possible stream names for the 
%                               "background" channel (e.g. {'x05A','405A'}).
%                               NOTE: Only one stream per file can be
%                               treated as background.
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
%       % Suppose you have 3 TDT blocks for which you want to extract data:
%       rawPaths       = ["C:\Data\Rat1\Block-1", ...
%                         "C:\Data\Rat1\Block-2", ...
%                         "C:\Data\Rat1\Block-3"];
%       extractedPaths = ["C:\Data\Rat1\Extracted\Block-1", ...
%                         "C:\Data\Rat1\Extracted\Block-2", ...
%                         "C:\Data\Rat1\Extracted\Block-3"];
%       sigNames       = {'x65A','465A','x465A'};
%       baqNames       = {'x05A','405A','x405A'};
%
%       extractTDTdata(rawPaths, extractedPaths, sigNames, baqNames, ...
%           'trim', 10, 'skipexisting', 0);
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
    for eachfile = 1:length(rawfolderpaths)
        % Prepare a new block data structure for this file
        blockdata = struct();
        blockdata.RawFolderPath = char(rawfolderpaths(eachfile));
        blockdata.ExtractedFolderPath = char(extractedfolderpaths(eachfile));

        % Construct output file name (where the extracted data will be saved)
        extractedMatFile = strcat(blockdata.ExtractedFolderPath, '_extracted.mat');
        
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
                streams = string(fieldnames(alldata.streams)); % Find all stream names for the current file
                epocs = string(fieldnames(alldata.epocs)); % Find all epoc names for the current file
                
                % Find the signal and background stream names for the current file
                currsigname = char(streams(ismember(streams,sigstreamnames))); % Extract the name of the signal field for the current file
                currbaqname = char(streams(ismember(streams,baqstreamnames))); % Extract the name of the background field for the current file

                % Display warning if no matching signal or background was found
                if isempty(currsigname)
                    warning('No valid SIGNAL stream found for file #%d.', eachfile);
                end
                if isempty(currbaqname)
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

                % Adjust epoc event times for the trimmed data based on epoc onset
                for eachepoc = 1:length(epocs)
                    thisepoc = char(epocs(eachepoc));
                    blockdata.(thisepoc) = round(((alldata.epocs.(thisepoc).onset-alldata.time_ranges(1)).*blockdata.fs)-(trimsamples));
                end

                % Add a record of which function + parameters were used
                blockdata.extractFunction  = mfilename;
                blockdata.params = params;
            
                % Save the extracted data struct as <ExtractedFolderPath>_extracted.mat
                save(extractedMatFile,'-struct', 'blockdata');

            catch ME % % If there's an error (e.g., path not found), log it and continue.
                fprintf('ERROR extracting file #%d. Verify raw and extract file paths.\n', eachfile);
                fprintf('   Raw path: %s\n', blockdata.RawFolderPath);
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
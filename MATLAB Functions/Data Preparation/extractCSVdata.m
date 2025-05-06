function [] = extractCSVdata(rawfolderpaths,extractedfolderpaths,sigstreamname,baqstreamname,varargin)
% EXTRACTCSVDATA  Extract and save fiber photometry data from generic csv
%                 file format. See PASTa user guide for details.
%
%   EXTRACTCSVDATA(RAWFOLDERPATHS, EXTRACTEDFOLDERPATHS, SIGSTREAMNAME, 
%   BAQSTREAMNAME, 'PARAM1', VAL1, 'PARAM2', VAL2, ...) reads csv files
%   containing raw fiber photometry data from the folders specified in 
%   RAWFOLDERPATHS, extracts signal and 
%   background streams from the filenames specified in SIGSTREAMNAME and 
%   BAQSTREAMNAME, trims the data streams, and saves the results to 
%   EXTRACTEDFOLDERPATHS in MATLAB .mat files.
%   The final extracted data structure is saved as:  <ExtractedFolderPath>_extracted.mat
%
% REQUIRED INPUTS:
%       RAWFOLDERPATHS        - String array of raw block folder paths.
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
%       SIGSTREAMNAME         - String; Name of csv files containing the "signal" channel.
%                               (e.g. 'sig'). NOTE: Only one stream per file can be
%                               treated as signal.
%
%       BAQSTREAMNAME         - String; Name of csv files containing the "background" channel.
%                               (e.g. 'baq'). NOTE: Only one stream per file can be
%                               treated as signal.
%
%  OPTIONAL INPUT NAME-VALUE PAIRS:
%       'loadepocs'   - Logical; Set to 1 to load event epoch files in the
%                       'Raw Data' session folders.
%
%       'epocsnames'  - Cell array of file names of csv files containing
%                       event epochs to be loaded into the data structure.
%
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
%       sigstreamname  = 'sig';
%       baqstreamname  = 'baq';
%
%       extractCSVdata(rawPaths, extractedPaths, sigstreamname, baqstreamname,
%           'loadepocs',1,'epocsnames',{'injt', 'strt'},'trim', 10, 'skipexisting', 0);
%
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
    addParameter(p, 'loadepocs', defaultparameters.loadepocs, @(x) any(x == [0 1])); % loadepocs: input must be 0 or 1
    addParameter(p, 'epocsnames',{}); % if loadepocs = 1, optional input epocsnames must include at least one epoch file name.

    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;

    % Main display and function inputs
    if params.skipexisting == 0
        disp('EXTRACTCSVDATA: all blocks will be extracted.')
    elseif params.skipexisting == 1
        disp('EXTRACTCSVDATA: pre-extracted blocks will be skipped.')
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
                % Load recording parameters
                recordingparams = readtable(fullfile(blockdata.RawFolderPath, 'recordingparams.csv'));
                recordingparamsNames = recordingparams.Properties.VariableNames; 
    
                for eachrecordingparam = 1:length(recordingparamsNames)
                    currparam = char(recordingparamsNames(eachrecordingparam));
                    blockdata.(currparam) = recordingparams.(currparam);
                end
        
                % Load signal and background streams
                sig = readmatrix(fullfile(blockdata.RawFolderPath, append(sigstreamname,'.csv')))';
                baq = readmatrix(fullfile(blockdata.RawFolderPath, append(baqstreamname,'.csv')))';

                % Trim the data at the beginning and end, as specified
                trimsamples = round(params.trim*blockdata.fs); % Find number of samples to trim from beginning and end of session
                blockdata.sig = sig((trimsamples+1):(end-trimsamples)); % Add signal to the data structure and trim
                blockdata.baq = baq((trimsamples+1):(end-trimsamples)); % Add signal to the data structure and trim

                % OPTIONAL: Load epocs and adjust epoch event indexes for the trimmed data based on epoc onset
                if params.loadepocs == 1
                    for eachepoc = 1:length(params.epocsnames)
                        currepocname = char(params.epocsnames(eachepoc));
                        currepoc = readmatrix(fullfile(blockdata.RawFolderPath, append(currepocname,'.csv')));
                        blockdata.(currepocname) = currepoc - trimsamples;
                    end
                end

                % Add a record of which function + parameters were used
                blockdata.params.(mfilename) = params;
                blockdata.params.(mfilename).extractFunction = mfilename;
            
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
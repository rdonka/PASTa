function [] = extractTDTdata(rawfolderpaths, extractedfolderpaths, sigstreamnames, baqstreamnames, varargin)
% EXTRACTTDTDATA    Extracts and saves raw fiber photometry data collected
%                   via Tucker Davis Technologies program Synapse, which
%                   saves raw data into a block structure.
%
% INPUTS:
%   RAWFOLDERPATHS:     String array; contains the full paths to the folder 
%                       location of the raw data blocks to be extracted. 
%                       The string array should contain one column with 
%                       each full path in a separate row. If using the
%                       loadKeys function, this can be created from the
%                       experiment key. 
%                       For example:
%                       rawfolderpaths = string({experimentkey.RawFolderPath})';
%
%   EXTRACTEDFOLDERPATHS: String array; contains the full paths to the folder 
%                       location in which to save the extracted MATLAB 
%                       structs for each block to be extracted. The string 
%                       array should contain one column with each full path 
%                       in a separate row. If using the loadKeys function, 
%                       this can be created from the experiment key. 
%                       For example:
%                       extractedfolderpaths = string({experimentkey.ExtractedFolderPath})';
%
%   SIGSTREAMNAMES:     Cell array; contains strings with the names of
%                       the streams to be treated as signal.
%                       NOTE: Only one stream per file can be treated as
%                       signal. If different files have different signal
%                       stream names, include all signal stream name
%                       variations in the cell array.
%                       For example: sigstreamnames = {'x65A', '465A', 'x465A'};
%
%   BAQSTREAMNAMES:     Cell array; contains strings with the names of
%                       the streams to be treated as background. 
%                       NOTE: Only one stream per file can be loaded to
%                       background. If different files have different 
%                       background stream names, include all background 
%                       stream name variations in the cell array.
%                       For example: baqstreamnames = {'x05A', '405A', 'x405A'};
%
% OPTIONAL INPUTS:
%   CLIP:               Numeric; Specified number of seconds to clip at the 
%                       beginning and end of the session.
%                       Default: 5
%
%   SKIPEXISTING:       Numeric; A binary variable. If set to 1, blocks 
%                       that have already been extracted will be skipped.
%                       If set to 0, pre-existing extracted blocks will be
%                       re-extracted.
%                       Default: 1
%
% OUTPUTS:
%       Saves the blockdata data structure to the ExtractedFolderPath with
%       the naming convention "ORIGINALFOLDERNAME_extracted.mat'
%
% Written by R M Donka, August 2024
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/
disp('BIN TRANSIENTS: Add bin variable to transient quantification table.')

    %% Prepare Settings
    % Import required and optional inputs into a structure
    inputs = struct(...
        'clip',[],...
        'skipexisting',[]);
    inputs = parseArgsLite(varargin,inputs);

    if isempty(inputs.clip) == true
        clip = 5;
    else
        clip = inputs.clip;
    end

    if isempty(inputs.skipexisting) == true
        skipexisting = 1;
    else
        skipexisting = inputs.skipexisting;
    end

    if skipexisting == 0
        disp('EXTRACTING TDT DATA: all blocks will be extracted.')
    elseif skipexisting == 1
        disp('EXTRACTING TDT DATA: pre-extracted blocks will be skipped.')
    end
       

    disp('INPUTS:') % Display all input values
    disp(inputs)

    %% Extract Raw Data
    for eachfile = 1:length(rawfolderpaths)
        blockdata = struct(); % Prepare block data structure
        blockdata.RawFolderPath = char(rawfolderpaths(eachfile)); % Pull out raw folder path
        blockdata.ExtractedFolderPath = char(extractedfolderpaths(eachfile)); % Pull out extracted folder path
        
        if skipexisting == 1 && isfile(strcat(blockdata.ExtractedFolderPath, '_extracted.mat')) % If skipexisting is set to 1, skip already extracted files
            fprintf('Skipping file number: %.f \n',eachfile) % Display which file is being skipped
            continue
        else % Extract block data
            try
                fprintf('Extracting file number: %.f \n',eachfile) % Display which file is loading
                alldata=TDTbin2mat(blockdata.RawFolderPath); % Load in all data using TDTbin2mat - see function for details.
        
                streams = string(fieldnames(alldata.streams)); % Find all stream names
                epocs = string(fieldnames(alldata.epocs)); % Find all epoc names
                
                currsigname = char(streams(ismember(streams,sigstreamnames))); % Extract the name of the signal field for the current file
                currbaqname = char(streams(ismember(streams,baqstreamnames))); % Extract the name of the background field for the current file
             
                disp(append('     Signal stream: ',currsigname)) % Display current file signal field name
                disp(append('     Background stream: ',currbaqname)) % Display current file background field name

                k = strfind(alldata.info.blockname,'-'); % Helper variable to pull out subject
                blockdata(1).Subject = alldata.info.blockname(1:k(1)-1); % Add subject ID to data
        
                blockdata.date = alldata.info.date; % Add date to data 
                blockdata.sessionduration = alldata.info.duration; % Add session duration to data
                blockdata.starttime = alldata.info.utcStartTime;
                blockdata.stoptime = alldata.info.utcStopTime;
        
                blockdata.fs = alldata.streams.(currsigname).fs; % Add sampling rate for each file to double check for consistency
                
                clipsamples = round(clip*blockdata.fs); % Find number of samples to clip from beginning and end of session
                blockdata.sig = double(alldata.streams.(currsigname).data((clipsamples+1):round(end-clipsamples))); % Add signal to the data structure and clip
                blockdata.baq = double(alldata.streams.(currbaqname).data((clipsamples+1):round(end-clipsamples))); % Add background to the data structure and clip

                for eachepoc = 1:length(epocs) % Adjust each epoc for clipping; Uses epoc onset
                    thisepoc = char(epocs(eachepoc));
                    blockdata.(thisepoc) = round(((alldata.epocs.(thisepoc).onset-alldata.time_ranges(1)).*blockdata.fs)-(clipsamples));
                end
            
                save(strcat(blockdata.ExtractedFolderPath, '_extracted.mat'),'-struct', 'blockdata'); % Save the extracted data structure to the extracted folder path location
            catch ME % If errors occur, usually the issue is that files cannot be located. Check paths. Loop continues to the next file.
                fprintf('ERROR: File %s\n', eachfile);
                fprintf('Message: %s\n', ME.message);
                continue
            end
        end
    end
end


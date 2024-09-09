function [] = extractTDTdata(rawfolderpaths, extractedfolderpaths, sigstreamnames, baqstreamnames, varargin)
% EXTRACTFPDATA     Extracts and saves raw fiber photometry data from TDT blocks.
%
% INPUTS:
%   RAWFOLDERPATHS:     A string array containing the paths to the folder 
%                       location of the raw data blocks to be extracted. 
%                       The string array should contain one column with 
%                       each full path in a separate row.
%
%   EXTRACTEDFOLDERPATHS: A string array containing the paths to the folder 
%                       location in which to save the extracted MatLab 
%                       structs for each block to be extracted. The string 
%                       array should contain one column with each full path 
%                       in a separate row.
%
%   SIGSTREAMNAMES:     A cell array containing strings with the names of
%                       the streams to be treated as signal. Note that
%                       only one stream per file can be treated as signal.
%                       If different files have different stream names,
%                       include all stream names in the cell array.
%
%   BAQSTREAMNAMES:     A cell array containing strings with the names of
%                       the streams to be treated as background. Note that
%                       only one stream per file can be treated as background.
%                       If different files have different stream names,
%                       include all stream names in the cell array.
%
% OPTIONAL INPUTS:
%
%   CLIP:               Specified number of seconds to clip at the beginning 
%                       and end of the session.
%
%   SKIPEXISTING:       A binary variable containing a 0 if pre-existing
%                       extracted blocks should be re-extracted or a 1 if
%                       pre-existing extracted blocks should be skipped.
%
% OUTPUTS:
%       Saves the blockdata data structure to the ExtractedFolderPath with
%       the naming convention "ORIGINALFOLDERNAME_extracted.mat'
%
% Written by R M Donka, August 2023
% Updated by R M Donka, July 2024
% Stored in RoitmanPhotometry GitHub repository, see Wiki for detailed notes.
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
        disp('EXTRACTING BLOCK DATA - all blocks will be extracted.')
    elseif skipexisting == 1
        disp('EXTRACTING BLOCK DATA - pre-extracted blocks will be skipped.')
    end
       
    for eachfile = 1:length(rawfolderpaths)
        blockdata = struct();
        blockdata.RawFolderPath = char(rawfolderpaths(eachfile));
        blockdata.ExtractedFolderPath = char(extractedfolderpaths(eachfile));
        
        if skipexisting == 1 && isfile(strcat(blockdata.ExtractedFolderPath, '_extracted.mat'))
            fprintf('Skipping file number: %.f \n',eachfile) % Display which file is being skipped
            continue
        else
            try % Use try catch in case any errors occur - usually this would be because a file can't be found
                fprintf('Extracting file number: %.f \n',eachfile) % Display which file is loading
                alldata=TDTbin2mat(blockdata.RawFolderPath); % Load in all data using TDTbin2mat - see function for details.
        
                streams = string(fieldnames(alldata.streams)); % Find all stream names
                epocs = string(fieldnames(alldata.epocs)); % Find all epoc names
                
                currsigname = char(streams(ismember(streams,sigstreamnames)));
                currbaqname = char(streams(ismember(streams,baqstreamnames)));
             
                disp(append('     Signal stream: ',currsigname))
                disp(append('     Background stream: ',currbaqname))

                k = strfind(alldata.info.blockname,'-'); % Helper variable to pull out subject
                blockdata(1).Subject = alldata.info.blockname(1:k(1)-1); % Add subject ID to data
        
                blockdata.date = alldata.info.date; % Add date to data 
                blockdata.sessionduration = alldata.info.duration; % Add session duration to data
                blockdata.starttime = alldata.info.utcStartTime;
                blockdata.stoptime = alldata.info.utcStopTime;
        
                blockdata.fs = alldata.streams.(currsigname).fs; % Add sampling rate for each file to double check for consistency
                
                clipsamples = round(clip*blockdata.fs);
                blockdata.sig = double(alldata.streams.(currsigname).data((clipsamples+1):round(end-clipsamples))); % Add 465 to data
                blockdata.baq = double(alldata.streams.(currbaqname).data((clipsamples+1):round(end-clipsamples))); % Add 405 to data

                for eachepoc = 1:length(epocs)
                    thisepoc = char(epocs(eachepoc));
                    blockdata.(thisepoc) = round(((alldata.epocs.(thisepoc).onset-alldata.time_ranges(1)).*blockdata.fs)-(clipsamples));
                end
            
                save(strcat(blockdata.ExtractedFolderPath, '_extracted.mat'),'-struct', 'blockdata'); 
            catch ME % Print an error message and go to the next file if loading fails
                fprintf('ERROR: File %s\n', eachfile);
                fprintf('Message: %s\n', ME.message);
                continue;
            end
        end
    end
end


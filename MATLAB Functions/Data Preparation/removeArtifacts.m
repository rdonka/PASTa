function [data] = removeArtifacts(data,whichstream,whichfs)
% FINDSESSIONTRANSIENTS_BLMEAN  Finds transients for the whole session. 
%                               Pre-transient baselines are set to the mean
%                               of the pre-transient window.
%
%                               NOTE: This sub-function is called by the
%                               main function FINDSESSIONTRANSIENTS. If
%                               called outside the main function, users
%                               must specify all input values manually.
%
% Copyright (C) 2025 Rachel Donka. Licensed under the GNU General Public License v3.
%
% INPUTS:
%       DATA:           This is a structure that contains at least the data
%                       stream for artifact removal.
%
%       WHICHSTREAM:    A variable containing a string with the name of the 
%                       field containing the stream for artifact removal. 
%                       For example, 'sigfilt'.
%                     
%       WHICHFS:        The name of the field containing the sampling rate
%                       of the streams (fs).
%
% OUTPUTS:
%       DATA:           The original data structure with WHICHSTREAM_AR,
%                       numartifcacts, and artifacts table added.
%                       added.  
%                       WHICHSTREAM_AR: Data stream with artifact periods
%                           replaced with NaNs.
%                       NUMARTIFACTS: Number of artifacts identified and
%                           replaced with NaNs in the stream.
%                       ARTIFACTS: Table containing the artifact peak loc
%                           and value, and removal period start and end
%                           indexes.
%
% Written by R M Donka, February 2025
% Stored in the PASTa GitHub Repository, see the user guide for additional
% documentation: https://rdonka.github.io/PASTa/

%% Prepare Settings
% Import required and optional inputs into a structure
    inputs = struct(...
        'whichstream',whichstream,...
        'whichfs',whichfs);

    outlierthresholdk = 3; % Outlier detection threshold; k = 3 is a common choice for outlier detection
    artifactremovalwindow = 0.3; % Seconds to replace with NaNs before and after artifact
    artifactampthresholdsd_max = 8; % Artifact SD detection threshold - max artifacts
    artifactampthresholdsd_min = 8; % Artifact SD detection threshold - min artifacts
    bucketsizeSecs = 30; % Bucket size (seconds) for session amplitude and mean calculations

    % Display settings
    disp("REMOVEARTIFACTS: Removing artifacts from subtracted signal.")

    disp('INPUTS:') % Display all input values
    disp(inputs)
    
    %% Find transients
    for eachfile = 1:length(data)
        try
            % Prep variables
            fs = data(eachfile).(whichfs);
            bucketsizeSamples = floor(bucketsizeSecs*fs);
            nbuckets = ceil(length(data(eachfile).(whichstream))/bucketsizeSamples);

            % Initialize arrays
            bucketvalues = [];
            overallvalues = [];
            allartifacts = [];
            
            % Calculate stream mean and amplitude for each bucket
            for eachbucket = 1:nbuckets
                if eachbucket == 1
                    bucketstart = 1;
                else
                    bucketstart = (eachbucket-1)*bucketsizeSamples;
                end
                if bucketstart+bucketsizeSamples-1 < length(data(eachfile).(whichstream))
                    bucketend = bucketstart+bucketsizeSamples-1;
                else
                    bucketend = length(data(eachfile).(whichstream));
                end
                bucketvalues.mean(eachbucket) = mean(data(eachfile).(whichstream)(bucketstart:bucketend));
                bucketvalues.amp(eachbucket) = max(data(eachfile).(whichstream)(bucketstart:bucketend)) - min(data(eachfile).(whichstream)(bucketstart:bucketend));
            end

            % Determine if any bucket amplitudes are outliers
            medVal = median(bucketvalues.amp); % Find median amplitude
            absDevs = abs(bucketvalues.amp - medVal); % Find deviations from median for all bucket amplitudes
            madVal = median(absDevs); % Find median deviation from median for all bucket amplitudes
            
            thresholdMAD = medVal + outlierthresholdk * madVal;  % Prepare MAD threshold
            outlierMaskMAD = bucketvalues.amp > thresholdMAD;  % Identify bins that exceed  outlier threshold

            if sum(outlierMaskMAD) > 0 % If any bin amplitudes are outliers, perform artifact detection
                binsinclude = outlierMaskMAD == 0; % Only use the non-outlier bins in the mean and amplitude calculations
                overallvalues.mean = mean(bucketvalues.mean(binsinclude)); % Find average stream mean
                overallvalues.amp = mean(bucketvalues.amp(binsinclude)); % Find average stream amplitude
                overallvalues.ampsd = std(bucketvalues.amp(binsinclude)); % Find stream amplitude standard deviation

                currartifactthreshold_max = (overallvalues.mean) + (overallvalues.amp/2) + (artifactampthresholdsd_max*overallvalues.ampsd); % Prepare max artifact threshold
                currartifactthreshold_min = -1*((overallvalues.mean) + (overallvalues.amp/2) + (artifactampthresholdsd_min*overallvalues.ampsd)); % Prepare min artifact threshold
    
                % Find artifacts
                allmaxlocs = find(islocalmax(data(eachfile).(whichstream))); % Find all maxes in stream
                allminlocs = find(islocalmin(data(eachfile).(whichstream))); % Find all mins in stream
    
                allmaxvals = data(eachfile).(whichstream)(allmaxlocs); % Find all max values
                allminvals = data(eachfile).(whichstream)(allminlocs); % Find all min values
    
                allartifactmaxlocs = allmaxlocs(allmaxvals>currartifactthreshold_max); % Find all maxes greater than artifact threshold
                allartifactminlocs = allminlocs(allminvals<currartifactthreshold_min); % Find all mins less than artifact threshold
    
                allartifactlocs = [allartifactmaxlocs allartifactminlocs]; % Create array of all max and min artifacts
                nanstream = data(eachfile).(whichstream); % Prepare stream in which to replace artifacts with NaNs
    
                if isempty(allartifactlocs) % If no artifacts found
                    data(eachfile).(append(whichstream,'_AR')) = nanstream; % Add NaN stream
                    data(eachfile).numartifacts = 0; % Add number of artifacts
                    data(eachfile).artifacts = allartifacts; % Add artifact variables
                else
                    disp(append('   WARNING File ',num2str(eachfile),': ',num2str(length(allartifactlocs)),' artifacts found.'))
                    allartifacts.peakloc = allartifactlocs; % Prepare all artifact peak locs
                    allartifacts.peakval = [data(eachfile).(whichstream)(allartifactlocs)]; % Prepare all artifact peak vals
                    allartifacts.startloc = allartifacts.peakloc-floor(fs*artifactremovalwindow); % Prepare all artifact removal start locs
                    allartifacts.endloc = allartifacts.peakloc+floor(fs*artifactremovalwindow); % Prepare all artifact removal end locs
    
                    for eachartifact = 1:length(allartifacts.peakloc) % Remove each artifact from stream and replace with NaNs
                        nanstream(allartifacts.startloc(eachartifact):allartifacts.endloc(eachartifact)) = NaN; 
                    end
    
                    data(eachfile).(append(whichstream,'_AR')) = nanstream; % Add NaN stream 
                    data(eachfile).numartifacts = length(allartifacts.peakloc); % Add number of artifacts
                    data(eachfile).artifacts = allartifacts; % Add artifact variables
                end
            end
        catch ME
            fprintf('   ERROR: %s\n', ME.message);
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

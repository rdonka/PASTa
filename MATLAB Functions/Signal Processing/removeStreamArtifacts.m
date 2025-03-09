function [artifactremovaldata] = removeStreamArtifacts(datastream,whichstream,fs,varargin)
% REMOVESTREAMARTIFACTS  Detects and removes artifacts from fiber photometry data streams.
%
%   ARTIFACTREMOVALDATA = REMOVESTREAMARTIFACTS(DATASTREAM, WHICHSTREAM, FS, 'PARAM1', VAL1, ...)
%   detects artifacts based on amplitude deviations and replaces them with NaNs or mean values.
%
% REQUIRED INPUTS:
%       DATASTREAM      - Numeric array; The signal data for artifact removal.
%
%       WHICHSTREAM     - String; Name of the stream being processed.
%
%       FS              - Numeric; Sampling rate of the data stream (Hz).
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%       'outlierthresholdk'     - Numeric, Integer; Multiplier for outlier 
%                                 threshold detection. Default: 3
%
%       'artifactremovalwindow' - Numeric; Seconds before and after an 
%                                 artifact to replace with NaNs. Default: 0.3
%
%       'artifactampthresh_max' - Numeric; Standard deviation threshold for 
%                                 high artifacts. Default: 8
%
%       'artifactampthresh_min' - Numeric; Standard deviation threshold for 
%                                 low artifacts. Default: 8
%
%       'bucketsizeSecs'        - Numeric; Bucket time window length 
%                                 (seconds) for amplitude and mean calculations.
%
% OUTPUTS:
%       ARTIFACTREMOVALDATA - Struct containing:
%           WHICHSTREAM_AR  - Data stream with artifacts replaced by NaNs.
%           WHICHSTREAM_AM  - Data stream with artifacts replaced by mean values.
%           NUMARTIFACTS    - Number of detected artifacts.
%           ARTIFACTS       - Table of artifact details (locations, values, start/end indices).
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
    addParameter(p, 'outlierthresholdk', defaultparameters.outlierthresholdk, @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'})); % Outlier detection threshold; must be a positive integer
    addParameter(p, 'artifactremovalwindow', defaultparameters.artifactremovalwindow, @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'})); % Artifact window; must be a positive integer
    addParameter(p, 'artifactampthreshold_max', defaultparameters.artifactampthreshold_max, @(x) validateattributes(x, {'numeric'}, {'positive'})); % Max artifact threshold (SD); must be positive
    addParameter(p, 'artifactampthreshold_min', defaultparameters.artifactampthreshold_min, @(x) validateattributes(x, {'numeric'}, {'positive'})); % Min artifact threshold (SD); must be positive
    addParameter(p, 'bucketsizeSecs', defaultparameters.bucketsizeSecs, @(x) validateattributes(x, {'numeric'}, {'positive'})); % Bucket size for amplitude calculations; must be positive
    parse(p, varargin{:});

    % Retrieve parsed inputs into params structure
    params = p.Results;
    
    %% Remove Artifacts
    % Initialize Variables
    bucketsizeSamples = floor(params.bucketsizeSecs*fs);
    nbuckets = ceil(length(datastream)/bucketsizeSamples);
    bucketvalues = struct('mean', [], 'amp', []);
    overallvalues = struct();
    allartifacts = struct();
    
    % Calculate Stream Statistics per Bucket
    for eachbucket = 1:nbuckets
        bucketstart = max(1, (eachbucket - 1) * bucketsizeSamples); % Start index of the current bucket
        bucketend = min(bucketstart + bucketsizeSamples - 1, length(datastream)); % End index of the current bucket
        
        bucketvalues.mean(eachbucket) = mean(datastream(bucketstart:bucketend)); % Compute mean signal value for this bucket
        bucketvalues.amp(eachbucket) = max(datastream(bucketstart:bucketend)) - min(datastream(bucketstart:bucketend)); % Compute signal amplitude for this bucket
    end

    % Identify Outlier Buckets Based on Amplitude
    medVal = median(bucketvalues.amp); % Compute median amplitude across buckets
    madVal = median(abs(bucketvalues.amp - medVal)); % Compute median absolute deviation (MAD)
    thresholdMAD = medVal + params.outlierthresholdk * madVal;  % Define outlier detection threshold
    outlierMaskMAD = bucketvalues.amp > thresholdMAD;  % Identify buckets exceeding the threshold

    if any(outlierMaskMAD) % If any bin amplitudes are outliers, perform artifact detection
        validBuckets = ~outlierMaskMAD; % Exclude outlier buckets from further calculations
        overallvalues.mean = mean(bucketvalues.mean(validBuckets)); % Compute overall mean from non-outlier buckets
        overallvalues.amp = mean(bucketvalues.amp(validBuckets)); % Compute overall amplitude from non-outlier buckets
        overallvalues.ampsd = std(bucketvalues.amp(validBuckets)); % Compute standard deviation of amplitude

        % Define artifact detection thresholds
        maxThreshold = overallvalues.mean + (overallvalues.amp / 2) + (params.artifactampthreshold_max * overallvalues.ampsd); % Upper threshold
        minThreshold = -1 * (overallvalues.mean + (overallvalues.amp / 2) + (params.artifactampthreshold_min * overallvalues.ampsd)); % Lower threshold
        
        % Detect artifacts by identifying local maxima and minima exceeding thresholds
        allmaxlocs = find(islocalmax(datastream)); % Find all local maxima in the data stream
        allminlocs = find(islocalmin(datastream)); % Find all local minima in the data stream
        
        % Extract locations of detected artifacts
        allartifactlocs = [allmaxlocs(datastream(allmaxlocs) > maxThreshold); allminlocs(datastream(allminlocs) < minThreshold)];
        nanstream = datastream; % Copy data stream to replace artifacts with NaNs
        meanstream = datastream; % Copy data stream to replace artifacts with mean values

        if isempty(allartifactlocs)
            % If no artifacts are detected, return the original data stream
            artifactremovaldata.(append(whichstream, 'AR')) = nanstream;
            artifactremovaldata.(append(whichstream, 'AM')) = meanstream;
            artifactremovaldata.numartifacts = 0;
            artifactremovaldata.artifacts = struct();
        else
            % Store artifact details
            allartifacts.peakloc = allartifactlocs;
            allartifacts.peakval = datastream(allartifactlocs);
            allartifacts.startloc = allartifactlocs - floor(fs * params.artifactremovalwindow);
            allartifacts.endloc = allartifactlocs + floor(fs * params.artifactremovalwindow);
            
            % Remove artifacts by replacing values with NaNs or mean values
            for eachartifact = 1:length(allartifacts.peakloc)
                nanstream(allartifacts.startloc(eachartifact):allartifacts.endloc(eachartifact)) = NaN; % Replace artifact region with NaNs
                meanstream(allartifacts.startloc(eachartifact):allartifacts.endloc(eachartifact)) = overallvalues.mean; % Replace artifact region with mean value
            end
            
            artifactremovaldata.(append(whichstream, 'AR')) = nanstream; % Store NaN-replaced stream
            artifactremovaldata.(append(whichstream, 'AM')) = meanstream; % Store mean-replaced stream
            artifactremovaldata.numartifacts = length(allartifacts.peakloc); % Store artifact count
            artifactremovaldata.artifacts = allartifacts; % Store artifact details
        end
    else
        % If no outlier buckets are found, return the original data stream
        artifactremovaldata.(append(whichstream, 'AR')) = datastream;
        artifactremovaldata.(append(whichstream, 'AM')) = datastream;
        artifactremovaldata.numartifacts = 0;
        artifactremovaldata.artifacts = struct();
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
